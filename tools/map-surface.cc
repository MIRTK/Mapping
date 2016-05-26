/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "mirtk/UniformSurfaceMapper.h"                   // Tutte (1964)
#include "mirtk/WeightedLeastSquaresSurfaceMapper.h"      // Kent et al. (1991), Floater (1997)
#include "mirtk/IntrinsicParameterizationSurfaceMapper.h" // Meyer et al. (2002)
#include "mirtk/MeanValueSurfaceMapper.h"                 // Floater (2003)

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  cout << "\n";
  cout << "usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "This tool computes a mapping for each point on the surface of a given input shape\n";
  cout << "embedded in 3D space. The output is a (piecewise linear) function which assigns each\n";
  cout << "point on the surface of the input shape one or more values. In case of non-closed surfaces,\n";
  cout << "the output map can interpolate any values given on the boundary of the surface at the\n";
  cout << "interior points of the surface. More common use cases are to compute a bijective mapping\n";
  cout << "from one geometric shape to another geometric shape with identical topology. The resulting\n";
  cout << "map is a parameterization of the surface of the input shape. Such parameterization can be\n";
  cout << "used for texturing, object morphing, and surface registration.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Point set delineating the boundary of the map domain or name\n";
  cout << "           of primitive shape such as: \"disk\", \"square\", or \"sphere\".\n";
  cout << "  output   File path of output map. A piecewise linear map is stored as VTK file.\n";
  cout << "           Other maps are stored in a custom binary format.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -barycentric   Use spring constants based on generalized barycentric coordiantes.\n";
  cout << "  -mean-value    Use spring constants based on mean value coordinates.\n";
  cout << "  -conformal     Conformal surface map or as-conformal-as-possible volumetric map.\n";
  cout << "  -harmonic      Harmonic volumetric map.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -p <n>                Exponent of harmonic energy term. When non-positive, solve for an\n";
  cout << "                        approximate harmonic surface map using a spring network. (default: 0)\n";
  cout << "  -name <string>        Name of point data array used as fixed point map.  (default: tcoords)\n";
  cout << "  -mask <string>        Name of point data array used as fixed point mask. (default: boundary)\n";
  cout << "  -max-iterations <n>   Maximum no. of linear solver iterations. (default: size of problem)\n";
  PrintCommonOptions(cout);
  cout << "\n";
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of implemented surface mapping methods
enum MapSurfaceMethod
{
  MAP_Uniform,               ///< Uniform edge weights, Tutte's planar graph mapping
  MAP_WeightedLeastSquares,  ///< Edge weights inverse proportional to edge length
  MAP_MeanValue,             ///< Floater's mean value convex map
  MAP_Intrinsic,             ///< Intrinsic parameterization with boundary constraints
  MAP_NaturalConformal,      ///< Meyer's natural conformal map
  MAP_PHarmonic,             ///< Joshi's p-harmonic map
  MAP_LeastSquaresConformal, ///< Levy's least squares conformal map
  MAP_Spectral,              ///< Spectral surface map w/o boundary constraints
  MAP_Spherical              ///< Spherical surface map w/o boundary constraints
};

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  const char      *values_name = nullptr;
  const char      *mask_name   = nullptr;
  MapSurfaceMethod method      = MAP_MeanValue;
  int              niters      = 0;

  int    pharmonic_exponent = 2;   // Exponent of p-harmonic energy
  int    wls_exponent       = 1;   // Weighted least squares exponent
  double intrinsic_lambda   = 0.0; // Conformal energy weight of intrinsic parameterization
  double intrinsic_mu       = 1.0; // Authalic  energy weight of intrinsic parameterization

  for (ALL_OPTIONS) {
    // Point data arrays
    if      (OPTION("-name")) values_name = ARGUMENT;
    else if (OPTION("-mask")) mask_name   = ARGUMENT;
    // Surface mapping method
    else if (OPTION("-uniform")) {
      method = MAP_Uniform;
    }
    else if (OPTION("-edge-length-weighted") || OPTION("-weighted-least-squares")) {
      method = MAP_WeightedLeastSquares;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(wls_exponent);
      else wls_exponent = 1;
    }
    else if (OPTION("-p-harmonic") || OPTION("-pharmonic")) {
      method = MAP_PHarmonic;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(pharmonic_exponent);
      else pharmonic_exponent = 2;
    }
    else if (OPTION("-intrinsic")) {
      method = MAP_Intrinsic;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(intrinsic_lambda);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(intrinsic_mu);
        else {
          intrinsic_lambda = clamp(intrinsic_lambda, .0, 1.0);
          intrinsic_mu     = 1.0 - intrinsic_lambda;
        }
      } else {
        intrinsic_lambda = 0.0;
        intrinsic_mu     = 1.0;
      }
    }
    else if (OPTION("-mean-value")) {
      method = MAP_MeanValue;
    }
    // Parameters of mapping method
    else if (OPTION("-max-iterations") || OPTION("-max-iter") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(niters);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  vtkSmartPointer<vtkPolyData>  mesh = ReadPolyData(input_name);
  vtkSmartPointer<vtkDataArray> values, mask;
  SharedPtr<Mapping>            map;

  vtkPointData * const pd = mesh->GetPointData();

    if (values_name) {
      values = pd->GetArray(values_name);
      if (values == nullptr) {
        FatalError("Input point set has no data array named " << values_name);
    } else {
      values = pd->GetTCoords();
      if (values == nullptr) {
        FatalError("Input point set has no TCOORDS, use -name option to specify name of map values array!");
      }
    }
  }

  if (mask_name) {
    mask = pd->GetArray(mask_name);
    if (mask == nullptr) {
      FatalError("Input point set has no data array named " << mask_name);
    }
  }

  switch (method) {
    case MAP_Uniform: {
      const char *msg = "Computing uniform surface map...";
      if (verbose) cout << msg, cout.flush();
      UniformSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Mesh(mesh);
      mapper.Input(values);
      mapper.Mask(mask);
      mapper.Run();
      map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_WeightedLeastSquares: {
      const char *msg = "Computing edge length weighted least squares surface map...";
      if (verbose) cout << msg, cout.flush();
      WeightedLeastSquaresSurfaceMapper mapper(wls_exponent);
      mapper.NumberOfIterations(niters);
      mapper.Mesh(mesh);
      mapper.Input(values);
      mapper.Mask(mask);
      mapper.Run();
      map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_MeanValue: {
      const char *msg = "Computing surface map using mean value coordinates...";
      if (verbose) cout << msg, cout.flush();
      MeanValueSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Mesh(mesh);
      mapper.Input(values);
      mapper.Mask(mask);
      mapper.Run();
      map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_Intrinsic: {
      const char *msg = "Computing surface map using intrinsic parameterization...";
      if (verbose) cout << msg, cout.flush();
      IntrinsicParameterizationSurfaceMapper mapper(intrinsic_lambda, intrinsic_mu);
      if (verbose) {
        cout << "\n  Conformal energy weight           = " << mapper.ConformalEnergyWeight();
        cout << "\n  Authalic  energy weight           = " << mapper.AuthalicEnergyWeight();
        cout.flush();
      }
      mapper.NumberOfIterations(niters);
      mapper.Mesh(mesh);
      mapper.Input(values);
      mapper.Mask(mask);
      mapper.Run();
      map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_PHarmonic: {
      FatalError("p-harmonic mapping using finite element method (FEM) not implemented");
      if (verbose) cout << "Computing p=" << pharmonic_exponent << " harmonic surface map...", cout.flush();
    } break;

    default: {
      FatalError("Selected mapping method not implemented");
    } break;
  }

  if (!map->Write(output_name)) {
    if (verbose) cout << " failed" << endl;
    FatalError("Failed to write surface map to " << output_name);
  }
  if (verbose) cout << " done" << endl;

  return 0;
}

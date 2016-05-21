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
#include "mirtk/LinearSpringSurfaceMapper.h"

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
  MAP_Barycentric, ///< Compute surface map using generalized barycentric coordinates
  MAP_MeanValue,   ///< Compute surface map using mean value coordinates
  MAP_Conformal,   ///< Compute conformal surface map
  MAP_Harmonic,    ///< Compute harmonic surface map
  MAP_PHarmonic,   ///< Comptue p-harmonic surface map
  MAP_Spectral,    ///< Compute surface map using spectral coordinates
  MAP_Spherical    ///< Compute spherical map of genus-0 surface
};

// -----------------------------------------------------------------------------
/// Compose planar surface map with inverse stereographic projection to sphere
vtkSmartPointer<vtkPolyData>
InverseStereographicProjection(vtkSmartPointer<vtkPolyData> surface)
{
  double bounds[6], extent[3], center[3], p[3], scale, radius, radius2;
  surface->GetBounds(bounds);
  extent[0] = bounds[1] - bounds[0];
  extent[1] = bounds[3] - bounds[2];
  extent[2] = bounds[5] - bounds[4];
  center[0] = bounds[0] + .5 * extent[0];
  center[1] = bounds[2] + .5 * extent[1];
  center[2] = bounds[4] + .5 * extent[2];
  radius    = .5 * max(max(extent[0], extent[1]), extent[2]);
  radius2   = radius * radius;
  vtkSmartPointer<vtkPoints> points = surface->GetPoints()->NewInstance();
  points->SetNumberOfPoints(surface->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    surface->GetPoint(ptId, p);
    p[0] -= center[0], p[1] -= center[1], p[2] -= center[2];
    scale = 2.0 * radius2 / (p[0]*p[0] + p[1]*p[1] + radius2);
    p[0] *= scale, p[1] *= scale, p[2] = (1.0 - scale) * radius;
    points->SetPoint(ptId, p);
  }
  vtkSmartPointer<vtkPolyData> output = surface->NewInstance();
  output->ShallowCopy(surface);
  output->SetPoints(points);
  return output;
}

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
  MapSurfaceMethod method      = MAP_Harmonic;
  int              niters      = 0;
  int              p           = 0; // exponent of p-harmonic energy functional
                                    // <=0: Use Eck's harmonic map approximation
                                    //  >0: Use A. Joshi's p-harmonic mapping method

  for (ALL_OPTIONS) {
    // Point data arrays
    if      (OPTION("-name")) values_name = ARGUMENT;
    else if (OPTION("-mask")) mask_name   = ARGUMENT;
    // Surface mapping method
    else if (OPTION("-barycentric"))  method = MAP_Barycentric;
    else if (OPTION("-mean-value"))   method = MAP_MeanValue;
    else if (OPTION("-conformal"))    method = MAP_Conformal;
    else if (OPTION("-harmonic"))     method = MAP_Harmonic;
    // Parameters of mapping method
    else if (OPTION("-p")) PARSE_ARGUMENT(p);
    else if (OPTION("-max-iterations") || OPTION("-max-iter") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(niters);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (method == MAP_Harmonic && p > 0) method = MAP_PHarmonic;

  vtkSmartPointer<vtkPolyData>  input = ReadPolyData(input_name);
  vtkSmartPointer<vtkDataArray> values, mask;
  SharedPtr<Mapping>            map;

  vtkPointData * const pd = input->GetPointData();

  if (values_name) {
    values = pd->GetArray(values_name);
    if (values == nullptr) {
      FatalError("Input point set has no data array named " << values_name);
    }
  } else {
    values = pd->GetTCoords();
    if (values == nullptr) {
      FatalError("Input point set has no TCOORDS, use -name option to specify name of map values array!");
    }
  }

  if (mask_name) {
    mask = input->GetPointData()->GetArray(mask_name);
    if (mask == nullptr) {
      FatalError("Input point set has no data array named " << mask_name);
    }
  }

  switch (method) {
    case MAP_Barycentric:
    case MAP_MeanValue:
    case MAP_Harmonic: {
      const char *what;
      LinearSpringSurfaceMapper::SpringConstantType type;
      if (method == MAP_Barycentric) {
        what = "barycentric coordinates";
        type = LinearSpringSurfaceMapper::BarycentricCoordinates;
      } else if (method == MAP_MeanValue) {
        what = "mean value coordinates";
        type = LinearSpringSurfaceMapper::MeanValueCoordinates;
      } else if (method == MAP_Harmonic) {
        what = "harmonic spring constants";
        type = LinearSpringSurfaceMapper::Harmonic;
      } else {
        what = "default spring constants";
        type = LinearSpringSurfaceMapper::Default;
      }
      if (verbose) cout << "Computing surface map using " << what << "...", cout.flush();
      typedef LinearSpringSurfaceMapper Mapper;
      Mapper mapper(type);
      mapper.NumberOfIterations(niters);
      mapper.Domain(input);
      mapper.Input(values);
      mapper.Mask(mask);
      mapper.Run();
      map = mapper.Output();
      if (verbose) cout << "Computing surface map using " << what << "...", cout.flush();
    } break;

    case MAP_PHarmonic: {
      FatalError("p-harmonic mapping using finite element method (FEM) not implemented");
      if (verbose) cout << "Computing p=" << p << " harmonic surface map...", cout.flush();
    } break;

    case MAP_Conformal: {
      FatalError("Conformal surface mapping not implemented");
    } break;

    case MAP_Spectral: {
      FatalError("Spectral mapping of a surface mesh not implemented");
    } break;

    case MAP_Spherical: {
      FatalError("Spherical mapping of a surface mesh not implemented");
    } break;
  }

  if (!map->Write(output_name)) {
    if (verbose) cout << " failed" << endl;
    FatalError("Failed to write surface map to " << output_name);
  }
  if (verbose) cout << " done" << endl;

  return 0;
}

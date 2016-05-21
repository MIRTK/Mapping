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

#include "mirtk/Vector.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/VtkMath.h"

#include "mirtk/BoundaryMapper.h"

#include "mirtk/BoundaryToDiskMapper.h"
#include "mirtk/BoundaryToSquareMapper.h"
#include "mirtk/BoundaryToPolygonMapper.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkUnsignedCharArray.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  cout << "\n";
  cout << "usage: " << name << " <domain> <codomain> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Compute a surface boundary map suitable for the computation of a\n";
  cout << "  surface map to a 2D primitive shape such as a disk, square, or polygon.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  domain     Tesselation of the surface of the input shape.\n";
  cout << "  codomain   File name of piecewise linear complex or name of primitive shape:\n";
  cout << "             - disk (default)\n";
  cout << "             - square\n";
  cout << "             - polygon\n";
  cout << "  output     Output point set with boundary values as point data.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -radius <float>\n";
  cout << "      Radius of primitive codomain shape. When non-positive, the radius is\n";
  cout << "      chosen such that the area matches the area of the input surface. (default: 0)\n";
  cout << "  -name <string>\n";
  cout << "      Name of boundary values point data array. (default: Map)\n";
  cout << "  -mask <string>\n";
  cout << "      Name of output boundary mask point data array. (default: none)\n";
  PrintCommonOptions(cout);
  cout << "\n";
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(3);

  const char *domain_name   = POSARG(1); // File name of input point set
  const char *codomain_name = POSARG(2); // Shape of codomain (boundary)
  const char *output_name   = POSARG(3); // File name of output point set
  const char *values_name   = "Map";     // Name of boundary values array
  const char *mask_name     = nullptr;   // Name of boundary mask array
  double      radius        = .0;        // Radius of primitive shape

  vtkSmartPointer<vtkIdList> selection = vtkSmartPointer<vtkIdList>::New();

  for (ALL_OPTIONS) {
    if      (OPTION("-name")) values_name = ARGUMENT;
    else if (OPTION("-mask")) mask_name = ARGUMENT;
    else if (OPTION("-codomain")) codomain_name = ARGUMENT;
    else if (OPTION("-radius")) PARSE_ARGUMENT(radius);
    else if (OPTION("-select")) {
      unsigned int ptId;
      do {
        PARSE_ARGUMENT(ptId);
        selection->InsertNextId(static_cast<vtkIdType>(ptId));
      } while (HAS_ARGUMENT);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  string codomain_shape = ToLower(codomain_name);
  if      (codomain_shape == "circle") codomain_shape = "disk";
  else if (codomain_shape == "poly")   codomain_shape = "polygon";
  else if (codomain_shape == "ngon")   codomain_shape = "polygon";
  else if (codomain_shape == "n-gon")  codomain_shape = "polygon";
  else if (codomain_shape == "ball")   codomain_shape = "sphere";

  if (codomain_shape == "sphere") {
    FatalError("A sphere has no surface boundary to be mapped!");
  }

  vtkSmartPointer<vtkPolyData> codomain;
  if (codomain_shape != "disk" && codomain_shape != "square" && codomain_shape != "polygon") {
    if (verbose) cout << "Reading codomain shape...", cout.flush();
    codomain = ReadPolyData(codomain_name);
    int nbounds, ncomps;
    double genus = Genus(codomain, nullptr, nullptr, nullptr, &nbounds, &ncomps, nullptr);
    if (genus != .0 || ncomps != 1 || nbounds != 1) {
      if (verbose) cout << " failed" << endl;
      FatalError("Surface map -codomain must be a single non-closed genus-0 surface!");
    }
    if (verbose) cout << " done" << endl;
    codomain_shape.clear();
  }

  if (verbose) cout << "Reading input surface...", cout.flush();
  vtkSmartPointer<vtkPointSet>  input    = ReadPointSet(domain_name);
  vtkSmartPointer<vtkPolyData>  domain   = DataSetSurface(input, true);
  vtkSmartPointer<vtkDataArray> orig_ids = domain->GetPointData()->GetArray("vtkOriginalPointIds");

  int nbounds, ncomps;
  double genus = Genus(domain, nullptr, nullptr, nullptr, &nbounds, &ncomps, nullptr);
  if (genus != .0 || ncomps != 1 || nbounds != 1) {
    FatalError("Input must be a single non-closed genus-0 surface!");
  }

  if (verbose) cout << " done\nComputing boundary map...", cout.flush();
  SharedPtr<BoundaryMapper> mapper;

  if (codomain_shape == "disk") {
    SharedPtr<BoundaryToDiskMapper> disk_mapper(new BoundaryToDiskMapper());
    disk_mapper->Radius(radius);
    mapper = disk_mapper;
  } else if (codomain_shape == "square") {
    SharedPtr<BoundaryToSquareMapper> square_mapper(new BoundaryToSquareMapper());
    square_mapper->SideLength(2.0 * radius);
    mapper = square_mapper;
  } else if (codomain_shape == "polygon") {
    SharedPtr<BoundaryToPolygonMapper> polygon_mapper(new BoundaryToPolygonMapper());
    polygon_mapper->Radius(radius);
    mapper = polygon_mapper;
  } else {
    FatalError("Mapping of surface boundary to arbitrary codomain polygon not implemented");
  }

  mapper->Surface(domain);
  mapper->Selection(selection);
  mapper->Run();

  vtkSmartPointer<vtkDataArray> mask;
  if (mask_name) {
    mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
    mask->SetName(mask_name);
    mask->SetNumberOfComponents(1);
    mask->SetNumberOfTuples(input->GetNumberOfPoints());
    mask->FillComponent(0, .0);
    input->GetPointData()->AddArray(mask);
  }

  vtkSmartPointer<vtkDataArray> values;
  vtkSmartPointer<vtkDataArray> bvalues = mapper->Values();
  values = vtkSmartPointer<vtkFloatArray>::New();
  values->SetName(values_name);
  values->SetNumberOfComponents(bvalues->GetNumberOfComponents());
  values->SetNumberOfTuples(input->GetNumberOfPoints());
  for (int j = 0; j < values->GetNumberOfComponents(); ++j) {
    values->FillComponent(j, .0);
  }
  for (vtkIdType i = 0, ptId; i < bvalues->GetNumberOfTuples(); ++i) {
    ptId = static_cast<vtkIdType>(mapper->BoundaryPointId(i));
    ptId = static_cast<vtkIdType>(orig_ids->GetComponent(ptId, 0));
    values->SetTuple(ptId, bvalues->GetTuple(i));
    if (mask) mask->SetComponent(ptId, 0, 1.0);
  }
  input->GetPointData()->SetTCoords(values);

  if (verbose) cout << " done" << endl;
  return WritePointSet(output_name, input) ? 0 : 1;
}
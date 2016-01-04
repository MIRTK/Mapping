/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#include <mirtkDiscreteMap.h>
#include <mirtkPointSetUtils.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkCellLocator.h>


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

const double DiscreteMap::_Tolerance2 = 1e-9;

// -----------------------------------------------------------------------------
void DiscreteMap::CopyAttributes(const DiscreteMap &other)
{
  _Input  = other._Input;
  _Values = other._Values;
}

// -----------------------------------------------------------------------------
DiscreteMap::DiscreteMap()
{
}

// -----------------------------------------------------------------------------
DiscreteMap::DiscreteMap(const DiscreteMap &other)
:
  VolumetricMap(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
DiscreteMap &DiscreteMap::operator =(const DiscreteMap &other)
{
  if (this != &other) {
    VolumetricMap::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void DiscreteMap::Initialize()
{
  // Check input mesh
  if (!_Input) {
    cerr << "DiscreteMap::Initialize: No input mesh set" << endl;
    exit(1);
  }
  if (_Input->GetNumberOfCells() == 0) {
    cerr << "DiscreteMap::Initialize: Input mesh has no cells" << endl;
    exit(1);
  }
  // Get point data array if not set
  if (!_Values) {
    if (_Input->GetPointData()->GetScalars()) {
      _Values = _Input->GetPointData()->GetScalars();
    } else if (_Input->GetPointData()->GetTCoords()) {
      _Values = _Input->GetPointData()->GetTCoords();
    } else if (_Input->GetPointData()->GetVectors()) {
      _Values = _Input->GetPointData()->GetVectors();
    } else {
      cerr << "DiscreteMap::Initialize: Input mesh has neither SCALARS, VECTORS, nor TCOORDS as point data!" << endl;
      exit(1);
    }
  }
  // Build cell locator
  _Locator = vtkSmartPointer<vtkCellLocator>::New();
  _Locator->SetDataSet(_Input);
  _Locator->BuildLocator();
}

// -----------------------------------------------------------------------------
VolumetricMap *DiscreteMap::NewCopy() const
{
  return new DiscreteMap(*this);
}

// -----------------------------------------------------------------------------
DiscreteMap::~DiscreteMap()
{
}

// =============================================================================
// Input domain
// =============================================================================

// -----------------------------------------------------------------------------
void DiscreteMap::BoundingBox(double &x1, double &y1, double &z1,
                                  double &x2, double &y2, double &z2) const
{
  if (_Input) {
    double bounds[6];
    _Input->GetBounds(bounds);
    x1 = bounds[0], x2 = bounds[1];
    y1 = bounds[2], y2 = bounds[3];
    z1 = bounds[4], z2 = bounds[5];
  } else {
    x1 = x2 = y1 = y2 = z1 = z2 = .0;
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
int DiscreteMap::NumberOfComponents() const
{
  return _Values->GetNumberOfComponents();
}

// -----------------------------------------------------------------------------
bool DiscreteMap::Evaluate(double *v, double x, double y, double z) const
{
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  double p[3] = {x, y, z}, pcoords[3];
  double *weight = new double[_Input->GetMaxCellSize()];
  if (_Locator->FindCell(p, _Tolerance2, cell, pcoords, weight) == -1) {
      for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
        v[j] = _OutsideValue;
      }
      delete[] weight;
      return false;
  }
  for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
    v[j] = .0;
  }
  vtkIdList * const ptIds = cell->GetPointIds();
  for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
    for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
      v[j] += weight[i] * _Values->GetComponent(ptIds->GetId(i), j);
    }
  }
  delete[] weight;
  return true;
}

// -----------------------------------------------------------------------------
double DiscreteMap::Evaluate(double x, double y, double z, int l) const
{
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  double p[3] = {x, y, z}, pcoords[3];
  double *weight = new double[_Input->GetMaxCellSize()];
  if (_Locator->FindCell(p, _Tolerance2, cell, pcoords, weight) == -1) {
      delete[] weight;
      return _OutsideValue;
  }
  double value = .0;
  vtkIdList * const ptIds = cell->GetPointIds();
  for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
    value += weight[i] * _Values->GetComponent(ptIds->GetId(i), l);
  }
  delete[] weight;
  return value;
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool DiscreteMap::Read(const char *fname)
{
  _Input = ReadPointSet(fname);
  if (_Input->GetNumberOfPoints() == 0) return false;
  if (_Input->GetPointData()->GetNumberOfArrays() == 1) {
    _Values = _Input->GetPointData()->GetArray(0);
  } else {
    _Values = NULL;
  }
  this->Initialize();
  return true;
}

// -----------------------------------------------------------------------------
bool DiscreteMap::Write(const char *fname) const
{
  vtkSmartPointer<vtkPointSet> output;
  output = vtkSmartPointer<vtkPointSet>::NewInstance(_Input);
  output->ShallowCopy(_Input);
  output->GetCellData ()->Initialize();
  output->GetPointData()->Initialize();
  if (_Values->GetNumberOfComponents() == 1) {
    output->GetPointData()->SetScalars(_Values);
  } else if (_Values->GetNumberOfComponents() == 3) {
    output->GetPointData()->SetTCoords(_Values);
  } else {
    output->GetPointData()->AddArray(_Values);
  }
  return WritePointSet(fname, output);
}


} // namespace mirtk

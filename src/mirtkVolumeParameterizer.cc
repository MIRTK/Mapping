/*
 * Medical Image Registration ToolKit (MMIRTK)
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

#include <mirtkVolumeParameterizer.h>

#include <mirtkMemory.h>
#include <mirtkPointSetUtils.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void VolumeParameterizer::CopyAttributes(const VolumeParameterizer &other)
{
  Delete(_OutputMap);
  _InputSet    = other._InputSet;
  _InputMap    = other._InputMap;
  _Boundary    = other._Boundary;
  _BoundaryMap = other._BoundaryMap;
  if (other._OutputMap) _OutputMap = other._OutputMap->NewCopy();
}

// -----------------------------------------------------------------------------
VolumeParameterizer::VolumeParameterizer()
:
  _OutputMap(NULL)
{
}

// -----------------------------------------------------------------------------
VolumeParameterizer::VolumeParameterizer(const VolumeParameterizer &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
VolumeParameterizer &VolumeParameterizer::operator =(const VolumeParameterizer &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
VolumeParameterizer::~VolumeParameterizer()
{
  Delete(_OutputMap);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void VolumeParameterizer::Run()
{
  this->Initialize();
  this->Parameterize();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void VolumeParameterizer::Initialize()
{
  // Free previous output map
  Delete(_OutputMap);

  // Check input
  if (!_InputSet) {
    cerr << "VolumeParameterizer::Initialize: Missing input point set" << endl;
    exit(1);
  }
  if (_InputSet->GetNumberOfCells() == 0) {
    cerr << "VolumeParameterizer::Initialize: Input has no cells" << endl;
    exit(1);
  }
  if (!_InputMap) {
    cerr << "VolumeParameterizer::Initialize: Missing input map" << endl;
    exit(1);
  }
  if (_InputMap->GetNumberOfTuples() != _InputSet->GetNumberOfPoints()) {
    cerr << "VolumeParameterizer::Initialize: Invalid input map" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void VolumeParameterizer::InitializeBoundary(vtkPointSet *input, vtkDataArray *map)
{
  vtkSmartPointer<vtkPointSet> volume;
  volume = vtkSmartPointer<vtkPointSet>::NewInstance(input);
  volume->ShallowCopy(input);
  volume->GetCellData ()->Initialize();
  volume->GetPointData()->Initialize();
  const int i  = volume->GetPointData()->AddArray(map);
  _Boundary    = DataSetSurface(volume, true);
  _BoundaryMap = _Boundary->GetPointData()->GetArray(i);
  _BoundaryMap->SetName("BoundaryMap");
}

// -----------------------------------------------------------------------------
void VolumeParameterizer::Finalize()
{
}


} // namespace mirtk

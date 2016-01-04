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

#include <mirtkTetrahedralVolumeParameterizer.h>

#include <mirtkVtk.h>
#include <mirtkPointSetUtils.h>
#include <mirtkDiscreteMap.h>

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
void TetrahedralVolumeParameterizer
::CopyAttributes(const TetrahedralVolumeParameterizer &other)
{
  _InputMask = other._InputMask;
  if (other._Volume && other._Coords && other._BoundaryMask) {
    _Coords = vtkSmartPointer<vtkDataArray>::NewInstance(other._Coords);
    _Coords->DeepCopy(other._Coords);
    _BoundaryMask = vtkSmartPointer<vtkDataArray>::NewInstance(other._BoundaryMask);
    _BoundaryMask->DeepCopy(other._BoundaryMask);
    _Volume = vtkSmartPointer<vtkPointSet>::NewInstance(other._Volume);
    _Volume->ShallowCopy(other._Volume);
    _Volume->GetPointData()->Initialize();
    _Volume->GetPointData()->AddArray(_Coords);
    _Volume->GetPointData()->AddArray(_BoundaryMask);
  } else {
    _Volume = NULL;
    _Coords = NULL;
    _BoundaryMask = NULL;
  }
  _NumberOfPoints         = other._NumberOfPoints;
  _NumberOfBoundaryPoints = other._NumberOfBoundaryPoints;
  _NumberOfInteriorPoints = other._NumberOfInteriorPoints;
}

// -----------------------------------------------------------------------------
TetrahedralVolumeParameterizer::TetrahedralVolumeParameterizer()
{
}

// -----------------------------------------------------------------------------
TetrahedralVolumeParameterizer
::TetrahedralVolumeParameterizer(const TetrahedralVolumeParameterizer &other)
:
  VolumeParameterizer(other),
  _NumberOfPoints(0),
  _NumberOfBoundaryPoints(0),
  _NumberOfInteriorPoints(0)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
TetrahedralVolumeParameterizer &TetrahedralVolumeParameterizer
::operator =(const TetrahedralVolumeParameterizer &other)
{
  if (this != &other) {
    VolumeParameterizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
TetrahedralVolumeParameterizer::~TetrahedralVolumeParameterizer()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void TetrahedralVolumeParameterizer::Initialize()
{
  // Initialize base class
  VolumeParameterizer::Initialize();

  // Tetrahedralize interior of input point set
  int map_index, mask_index = -1;
  vtkSmartPointer<vtkPointSet> input;
  input = vtkSmartPointer<vtkPointSet>::NewInstance(_InputSet);
  input->ShallowCopy(_InputSet);
  input->GetCellData ()->Initialize();
  input->GetPointData()->Initialize();
  map_index = input->GetPointData()->AddArray(_InputMap);
  if (_InputMask) mask_index = input->GetPointData()->AddArray(_InputMask);
  _Volume = Tetrahedralize(input);
  _Coords = _Volume->GetPointData()->GetArray(map_index);
  _Coords->SetName("VolumetricMap");

  // Extract surface of volume mesh
  this->InitializeBoundary(_Volume, _Coords);

  // Initialize boundary mask
  if (mask_index != -1) {
    _BoundaryMask = _Volume->GetPointData()->GetArray(mask_index);
  } else {
    _BoundaryMask = NewVTKDataArray(VTK_UNSIGNED_CHAR);
    _BoundaryMask->SetNumberOfComponents(1);
    _BoundaryMask->SetNumberOfTuples(_NumberOfPoints);
    _BoundaryMask->FillComponent(0, .0);
  }
  _BoundaryMask->SetName("BoundaryMask");
  vtkDataArray *origPtIds = _Boundary->GetPointData()->GetArray("vtkOriginalPointIds");
  for (vtkIdType ptId = 0, origPtId; ptId < _Boundary->GetNumberOfPoints(); ++ptId) {
    origPtId = static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0));
    _BoundaryMask->SetComponent(origPtId, 0, 1.0);
  }
  _NumberOfPoints         = static_cast<int>(_Volume->GetNumberOfPoints());
  _NumberOfBoundaryPoints = 0;
  for (vtkIdType ptId = 0; ptId < _BoundaryMask->GetNumberOfTuples(); ++ptId) {
    if (IsBoundaryPoint(ptId)) ++_NumberOfBoundaryPoints;
  }
  _NumberOfInteriorPoints = _NumberOfPoints - _NumberOfBoundaryPoints;
}

// -----------------------------------------------------------------------------
void TetrahedralVolumeParameterizer::Finalize()
{
  // Create volumetric output map
  DiscreteMap *map = new DiscreteMap();
  map->Input (_Volume);
  map->Values(_Coords);
  _OutputMap = map;

  // Finalize base class
  VolumeParameterizer::Finalize();
}


} // namespace mirtk

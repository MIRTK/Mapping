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

#include "mirtk/SurfaceMapper.h"

#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkUnsignedCharArray.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceMapper::CopyAttributes(const SurfaceMapper &other)
{
  _Domain      = other._Domain;
  _Surface     = other._Surface;
  _Input       = other._Input;
  _Mask        = other._Mask;
  _Fixed       = other._Fixed;
  _FixedPoints = other._FixedPoints;
  _FreePoints  = other._FreePoints;
  _PointIndex  = other._PointIndex;

  if (other._Values) {
    _Values = other._Values->NewInstance();
    _Values->DeepCopy(other._Values);
  } else {
    _Values = nullptr;
  }

  if (other._Output) {
    _Output = SharedPtr<Mapping>(other._Output->NewCopy());
  } else {
    _Output = nullptr;
  }
}

// -----------------------------------------------------------------------------
SurfaceMapper::SurfaceMapper()
{
}

// -----------------------------------------------------------------------------
SurfaceMapper::SurfaceMapper(const SurfaceMapper &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceMapper &SurfaceMapper::operator =(const SurfaceMapper &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceMapper::~SurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceMapper::Run()
{
  this->Initialize();
  this->Solve();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void SurfaceMapper::Initialize()
{
  // Free previous output map
  _Output = nullptr;

  // Check input
  if (!_Domain) {
    cerr << this->NameOfType() << "::Initialize: Missing input surface" << endl;
    exit(1);
  }
  if (_Domain->GetNumberOfPolys() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input point set must be a surface mesh" << endl;
    exit(1);
  }
  if (!_Input) {
    cerr << this->NameOfType() << "::Initialize: Missing boundary conditions" << endl;
    exit(1);
  }
  if (_Input && _Input->GetNumberOfTuples() != _Domain->GetNumberOfPoints()) {
    cerr << this->NameOfType() << "::Initialize: Invalid input map values array" << endl;
    exit(1);
  }
  if (_Fixed && _Fixed->GetNumberOfTuples() != _Domain->GetNumberOfPoints()) {
    cerr << this->NameOfType() << "::Initialize: Invalid input mask" << endl;
    exit(1);
  }

  // Make copy of input values array
  _Values = _Input->NewInstance();
  _Values->SetName(_Input->GetName());
  _Values->DeepCopy(_Input);

  // Internal surface mesh with point data arrays
  _Surface = _Domain->NewInstance();
  _Surface->ShallowCopy(_Domain);
  _Surface->SetLines(nullptr);
  _Surface->SetVerts(nullptr);
  _Surface->GetCellData()->Initialize();
  _Surface->GetPointData()->Initialize();
  _Surface->GetPointData()->SetTCoords(_Values);
  _Surface->GetPointData()->SetScalars(_Mask);

  // Remesh surface if necessary
  if (this->Remesh()) {
    _Values = _Surface->GetPointData()->GetTCoords();
  }

  // Build links
  _Surface->BuildLinks();

  // Initialize boundary mask
  _Fixed = _Surface->GetPointData()->GetScalars();
  if (!_Fixed) {
    _Fixed = vtkSmartPointer<vtkUnsignedCharArray>::New();
    _Fixed->SetNumberOfComponents(1);
    _Fixed->SetNumberOfTuples(_Surface->GetNumberOfPoints());
    _Fixed->FillComponent(0, .0);
    UnorderedSet<int> ptIds = BoundaryPoints(_Surface);
    for (auto it = ptIds.begin(); it != ptIds.end(); ++it) {
      _Fixed->SetComponent(static_cast<vtkIdType>(*it), 0, 1.0);
    }
  }

  // Set disjoint sets of IDs of free and fixed points
  _FixedPoints = vtkSmartPointer<vtkIdList>::New();
  _FixedPoints->Allocate(_Surface->GetNumberOfPoints());
  _FreePoints = vtkSmartPointer<vtkIdList>::New();
  _FreePoints->Allocate(_Surface->GetNumberOfPoints());
  _PointIndex = vtkSmartPointer<vtkIdList>::New();
  _PointIndex->SetNumberOfIds(_Surface->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < _Surface->GetNumberOfPoints(); ++ptId) {
    if (IsFixedPoint(ptId)) {
      _PointIndex->SetId(ptId, -(_FixedPoints->InsertNextId(ptId) + 1));
    } else {
      _PointIndex->SetId(ptId, _FreePoints->InsertNextId(ptId));
    }
  }
  _FixedPoints->Squeeze();
  _FreePoints ->Squeeze();
}

// -----------------------------------------------------------------------------
bool SurfaceMapper::Remesh()
{
  // Override in subclass if surface mesh must meet certain topological requirements
  return false;
}

// -----------------------------------------------------------------------------
void SurfaceMapper::Finalize()
{
  if (_Output == nullptr) {
    SharedPtr<PiecewiseLinearMap> map = NewShared<PiecewiseLinearMap>();
    map->Domain(_Surface);
    map->Values(_Values);
    _Output = map;
  }
}


} // namespace mirtk

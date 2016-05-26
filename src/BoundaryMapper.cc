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

#include "mirtk/BoundaryMapper.h"

#include "mirtk/Algorithm.h"
#include "mirtk/Memory.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ChordLengthBoundaryParameterizer.h"

#include "vtkFloatArray.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryMapper::CopyAttributes(const BoundaryMapper &other)
{
  _Surface   = other._Surface;
  _EdgeTable = other._EdgeTable;

  _BoundarySegments   = other._BoundarySegments;
  _BoundaryPointIds   = other._BoundaryPointIds;
  _BoundaryPointIndex = other._BoundaryPointIndex;

  _Parameterizer = SharedPtr<BoundaryParameterizer>(other._Parameterizer->NewCopy());

  if (other._Values) {
    _Values = other._Values->NewInstance();
    _Values->DeepCopy(other._Values);
  } else {
    _Values = nullptr;
  }
}

// -----------------------------------------------------------------------------
BoundaryMapper::BoundaryMapper()
:
  _Parameterizer(NewShared<ChordLengthBoundaryParameterizer>())
{
}

// -----------------------------------------------------------------------------
BoundaryMapper::BoundaryMapper(const BoundaryMapper &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryMapper &BoundaryMapper::operator =(const BoundaryMapper &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryMapper::~BoundaryMapper()
{
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
int BoundaryMapper::NumberOfComponents() const
{
  return 2;
}

// -----------------------------------------------------------------------------
int BoundaryMapper::NumberOfBoundarySegments() const
{
  if (_BoundarySegments) {
    return static_cast<int>(_BoundarySegments->GetNumberOfCells());
  } else if (_Surface) {
    return mirtk::NumberOfBoundarySegments(_Surface, _EdgeTable.get());
  } else {
    return 0;
  }
}

// -----------------------------------------------------------------------------
void BoundaryMapper::BoundaryPointIndices(int n, Array<int> &i) const
{
  vtkIdType npts, *pts;
  _BoundarySegments->GetCell(n, npts, pts);
  i.resize(npts);
  for (vtkIdType j = 0; j < npts; ++j) {
    i[j] = BoundaryPointIndex(pts[j]);
  }
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryMapper::Run()
{
  this->Initialize();
  for (int n = 0; n < NumberOfBoundarySegments(); ++n) {
    this->MapBoundary(n); 
  }
  this->Finalize();
}

// -----------------------------------------------------------------------------
void BoundaryMapper::Initialize()
{
  // Extract boundary segments
  if (!_EdgeTable) _EdgeTable.reset(new mirtk::EdgeTable(_Surface));
  _BoundarySegments = BoundarySegments(_Surface, _EdgeTable.get());
  if (NumberOfBoundarySegments() == 0) {
    cerr << this->NameOfType() << "::Initialize: Surface is closed and has no boundary segments to map!" << endl;
    exit(1);
  }

  // Obtain set of boundary point IDs
  _BoundaryPointIndex = vtkSmartPointer<vtkIdList>::New();
  _BoundaryPointIndex->SetNumberOfIds(_Surface->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < _BoundaryPointIndex->GetNumberOfIds(); ++ptId) {
    _BoundaryPointIndex->SetId(ptId, -1);
  }
  _BoundaryPointIds = vtkSmartPointer<vtkIdList>::New();
  _BoundaryPointIds->Allocate(_BoundarySegments->GetNumberOfConnectivityEntries());
  for (vtkIdType i = 0, npts, *pts; i < _BoundarySegments->GetNumberOfCells(); ++i) {
    _BoundarySegments->GetCell(i, npts, pts);
    for (vtkIdType j = 0; j < npts; ++j) {
      _BoundaryPointIndex->SetId(pts[j], _BoundaryPointIds->InsertNextId(pts[j]));
    }
  }
  _BoundaryPointIds->Squeeze();

  // Allocate boundary values array
  const int       dim = this->NumberOfComponents();
  const vtkIdType num = this->NumberOfBoundaryPoints();
  _Values = vtkSmartPointer<vtkFloatArray>::New();
  _Values->SetName("Map");
  if (_Values->GetNumberOfTuples()     != num ||
      _Values->GetNumberOfComponents() != dim) {
    _Values->SetNumberOfComponents(dim);
    _Values->SetNumberOfTuples(num);
  }
  for (int j = 0; j < dim; ++j) {
    _Values->FillComponent(j, .0);
  }
}

// -----------------------------------------------------------------------------
void BoundaryMapper::MapBoundary(int n)
{
  // Parameterize boundary segment
  _Parameterizer->Boundary(BoundaryPointIds(n));
  _Parameterizer->Selection(_Selection);
  _Parameterizer->Points(_Surface->GetPoints());
  _Parameterizer->Run();

  // Sort boundary points by increasing parameter value
  const int        npoints = NumberOfBoundaryPoints(n);
  const Array<int> ptIdx   = BoundaryPointIndices(n);

  Array<int> indices(npoints);
  for (int j = 0; j < npoints; ++j) indices[j] = j;
  SortIndicesOfArray<double> predicate(_Parameterizer->Values());
  sort(indices.begin(), indices.end(), predicate);

  Array<int>    i(npoints);
  Array<double> t(npoints);
  Array<int>    selection;
  selection.reserve(_Parameterizer->NumberOfSelectedPoints());
  for (int j = 0; j < npoints; ++j) {
    const auto &idx = indices[j];
    i[j] = ptIdx[idx];
    t[j] = _Parameterizer->Values()[idx];
    if (_Parameterizer->IsSelected(idx)) {
      selection.push_back(j);
    }
  }

  // Assign map values to points of boundary segment
  this->MapBoundarySegment(n, i, t, selection);
}

// -----------------------------------------------------------------------------
void BoundaryMapper::Finalize()
{
}


} // namespace mirtk

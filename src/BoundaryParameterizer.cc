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

#include "mirtk/BoundaryParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryParameterizer::CopyAttributes(const BoundaryParameterizer &other)
{
  _Boundary          = other._Boundary;
  _Points            = other._Points;
  _Selection         = other._Selection;
  _BoundarySelection = other._BoundarySelection;
  _Values            = other._Values;
}

// -----------------------------------------------------------------------------
BoundaryParameterizer::BoundaryParameterizer()
:
  _BoundarySelection(vtkSmartPointer<vtkIdList>::New())
{
}

// -----------------------------------------------------------------------------
BoundaryParameterizer::BoundaryParameterizer(const BoundaryParameterizer &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryParameterizer &BoundaryParameterizer::operator =(const BoundaryParameterizer &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryParameterizer::~BoundaryParameterizer()
{
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
Vector BoundaryParameterizer::BoundaryEdgeLengths() const
{
  const int npoints = NumberOfBoundaryPoints();
  Vector l(npoints);
  double p1[3], p2[3];
  BoundaryPoint(0, p1);
  for (int i = 1; i < npoints; ++i) {
    BoundaryPoint(i, p2);
    l(i-1) = Distance(p1, p2);
    p1[0] = p2[0], p1[1] = p2[1], p1[2] = p2[2];
  }
  BoundaryPoint(0, p2);
  l(npoints-1) = Distance(p1, p2);
  return l;
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryParameterizer::Run()
{
  this->Initialize();
  this->Parameterize();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void BoundaryParameterizer::Initialize()
{
  // Check boundary segment input
  if (!_Boundary) {
    cerr << this->NameOfType() << "::Initialize: No boundary points given!" << endl;
    exit(1);
  }
  if (!_Points) {
    cerr << this->NameOfType() << "::Initialize: No point set given!" << endl;
    exit(1);
  }
  const int npoints = NumberOfBoundaryPoints();
  const int maxPtId = static_cast<int>(_Points->GetNumberOfPoints()) - 1;
  if (npoints < 2) {
    cerr << this->NameOfType() << "::Initialize: Boundary segment must have at least two points!" << endl;
    exit(1);
  }
  for (int i = 0, ptId; i < npoints; ++i) {
    ptId = BoundaryPointId(i);
    if (ptId < 0 || ptId > maxPtId) {
      cerr << this->NameOfType() << "::Initialize: Invalid point index: " << ptId << endl;
      exit(1);
    }
  }
  // Ensure uniqueness of selected points and discard those that do not belong to this segment
  if (_Selection) {
    _BoundarySelection->Allocate(_Selection->GetNumberOfIds());
    for (vtkIdType i = 0, ptId; i < _Selection->GetNumberOfIds(); ++i) {
      ptId = _Selection->GetId(i);
      if (IsBoundaryPoint(ptId)) {
        _BoundarySelection->InsertUniqueId(_Selection->GetId(i));
      }
    }
    _BoundarySelection->Squeeze();
  } else {
    _BoundarySelection->Reset();
  }
  // Allocate parameter values
  _Values.resize(npoints);
}

// -----------------------------------------------------------------------------
void BoundaryParameterizer::Finalize()
{
  for (auto it = _Values.begin(); it != _Values.end(); ++it) {
    if (*it < .0 || *it >= 1.0) {
      cerr << this->NameOfType() << "::Finalize: Boundary curve parameter value must be in [0, 1), but t=" << *it << endl;
      exit(1);
    }
  }
}


} // namespace mirtk

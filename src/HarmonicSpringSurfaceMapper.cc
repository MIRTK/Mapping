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

#include "mirtk/HarmonicSpringSurfaceMapper.h"

#include "mirtk/PointSetUtils.h"
#include "mirtk/VtkMath.h"

#include "vtkIdList.h"
#include "vtkTriangle.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void HarmonicSpringSurfaceMapper::CopyAttributes(const HarmonicSpringSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
HarmonicSpringSurfaceMapper::HarmonicSpringSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
HarmonicSpringSurfaceMapper::HarmonicSpringSurfaceMapper(const HarmonicSpringSurfaceMapper &other)
:
  SymmetricLinearSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
HarmonicSpringSurfaceMapper &HarmonicSpringSurfaceMapper
::operator =(const HarmonicSpringSurfaceMapper &other)
{
  if (this != &other) {
    SymmetricLinearSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
HarmonicSpringSurfaceMapper::~HarmonicSpringSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
bool HarmonicSpringSurfaceMapper::Remesh()
{
  if (!IsTriangularMesh(_Surface)) {
    _Surface = Triangulate(_Surface);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
double HarmonicSpringSurfaceMapper::Weight(int i, int j) const
{
  vtkIdType npts, *pts;
  vtkIdType ptId1 = static_cast<vtkIdType>(i);
  vtkIdType ptId2 = static_cast<vtkIdType>(j);
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  double p1[3], p2[3], p3[3], l12, l13, l23, w = .0;
  _Surface->GetPoint(ptId1, p1);
  _Surface->GetPoint(ptId2, p2);
  _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
  l12 = vtkMath::Distance2BetweenPoints(p1, p2);
  for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
    _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
    while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
    _Surface->GetPoint(pts[0], p3);
    l13 = vtkMath::Distance2BetweenPoints(p1, p3);
    l23 = vtkMath::Distance2BetweenPoints(p2, p3);
    w += (l13 + l23 - l12) / vtkTriangle::TriangleArea(p1, p2, p3);
  }
  return w;
}


} // namespace mirtk

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

#include "mirtk/IntrinsicParameterizationSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"
#include "mirtk/PointSetUtils.h"

#include "vtkIdList.h"

namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void IntrinsicParameterizationSurfaceMapper
::CopyAttributes(const IntrinsicParameterizationSurfaceMapper &other)
{
  _ConformalEnergyWeight = other._ConformalEnergyWeight;
  _AuthalicEnergyWeight  = other._AuthalicEnergyWeight;
}

// -----------------------------------------------------------------------------
IntrinsicParameterizationSurfaceMapper
::IntrinsicParameterizationSurfaceMapper(double lambda)
:
  _ConformalEnergyWeight(clamp(lambda, .0, 1.0)),
  _AuthalicEnergyWeight(1.0 - _ConformalEnergyWeight)
{
}

// -----------------------------------------------------------------------------
IntrinsicParameterizationSurfaceMapper
::IntrinsicParameterizationSurfaceMapper(double lambda, double mu)
:
  _ConformalEnergyWeight(lambda),
  _AuthalicEnergyWeight(mu)
{
}

// -----------------------------------------------------------------------------
IntrinsicParameterizationSurfaceMapper
::IntrinsicParameterizationSurfaceMapper(
  const IntrinsicParameterizationSurfaceMapper &other
) :
  NonSymmetricLinearSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
IntrinsicParameterizationSurfaceMapper &
IntrinsicParameterizationSurfaceMapper
::operator =(const IntrinsicParameterizationSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricLinearSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
IntrinsicParameterizationSurfaceMapper
::~IntrinsicParameterizationSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
bool IntrinsicParameterizationSurfaceMapper::Remesh()
{
  if (!IsTriangularMesh(_Surface)) {
    _Surface = Triangulate(_Surface);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
double IntrinsicParameterizationSurfaceMapper::Weight(int i, int j) const
{
  vtkIdType npts, *pts;
  double p1[3], p2[3], p3[3], w = .0;
  vtkIdType ptId1 = static_cast<vtkIdType>(i);
  vtkIdType ptId2 = static_cast<vtkIdType>(j);
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  // Get edge points and adjacent cells
  _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
  _Surface->GetPoint(ptId1, p1);
  _Surface->GetPoint(ptId2, p2);
  // Dirichlet energy weights
  if (_ConformalEnergyWeight != .0) {
    for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
      _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
      while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
      _Surface->GetPoint(pts[0], p3);
      w += _ConformalEnergyWeight * Cotangent(p1, p3, p2);
    }
  }
  // Chi energy weights
  if (_AuthalicEnergyWeight != .0) {
    const double norm = _AuthalicEnergyWeight / vtkMath::Distance2BetweenPoints(p1, p2);
    for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
      _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
      while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
      _Surface->GetPoint(pts[0], p3);
      w += norm * Cotangent(p1, p2, p3);
    }
  }
  return w;
}


} // namespace mirtk

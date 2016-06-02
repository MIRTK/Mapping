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

#include "mirtk/MeanValueSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkPolyData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeanValueSurfaceMapper::CopyAttributes(const MeanValueSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::MeanValueSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::MeanValueSurfaceMapper(const MeanValueSurfaceMapper &other)
:
  NonSymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper &MeanValueSurfaceMapper::operator =(const MeanValueSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::~MeanValueSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double MeanValueSurfaceMapper::Weight(int i, int j) const
{
  double p1[3], p2[3], p3[3], e1[3], e2[3], l12, w = .0;

  vtkIdType npts, *pts;
  vtkIdType ptId1 = static_cast<vtkIdType>(i);
  vtkIdType ptId2 = static_cast<vtkIdType>(j);
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

  _Surface->GetPoint(ptId1, p1);
  _Surface->GetPoint(ptId2, p2);
  vtkMath::Subtract(p2, p1, e1);
  l12 = vtkMath::Normalize(e1);

  _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
  for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
    _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
    if (npts != 3) {
      cerr << this->NameOfType() << "::Weight: Surface mesh cells must be triangles" << endl;
      exit(1);
    }
    while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
    _Surface->GetPoint(pts[0], p3);
    vtkMath::Subtract(p3, p1, e2);
    vtkMath::Normalize(e2);
    w += tan(.5 * acos(vtkMath::Dot(e1, e2))) / l12;
  }

  return w;
}


} // namespace mirtk

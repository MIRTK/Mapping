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

#include "mirtk/LinearSpringSurfaceMapper.h"

#include "mirtk/Assert.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/VtkMath.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkTriangle.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void LinearSpringSurfaceMapper::CopyAttributes(const LinearSpringSurfaceMapper &other)
{
  _SpringConstant     = other._SpringConstant;
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
}

// -----------------------------------------------------------------------------
LinearSpringSurfaceMapper::LinearSpringSurfaceMapper(SpringConstantType type)
:
  _SpringConstant(type),
  _NumberOfIterations(0),
  _Tolerance(.0)
{
}

// -----------------------------------------------------------------------------
LinearSpringSurfaceMapper::LinearSpringSurfaceMapper(const LinearSpringSurfaceMapper &other)
:
  SurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LinearSpringSurfaceMapper &LinearSpringSurfaceMapper
::operator =(const LinearSpringSurfaceMapper &other)
{
  if (this != &other) {
    SurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LinearSpringSurfaceMapper::~LinearSpringSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
bool LinearSpringSurfaceMapper::Symmetric() const
{
  return _SpringConstant != BarycentricCoordinates &&
         _SpringConstant != MeanValueCoordinates;
}

// -----------------------------------------------------------------------------
// See Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
inline double Cotangent(double a[3], double b[3], double c[3])
{
  double ba[3], bc[3], n[3];
  vtkMath::Subtract(a, b, ba);
  vtkMath::Subtract(c, b, bc);
  vtkMath::Cross(ba, bc, n);
  return vtkMath::Dot(ba, bc) / vtkMath::Norm(n);
}

// -----------------------------------------------------------------------------
double LinearSpringSurfaceMapper::K(int i, int j) const
{
  double k_ij = .0;
  switch (_SpringConstant)
  {
    // Tutte (1963). How to draw a graph.
    case Uniform: {
      k_ij = 1.0;
    } break;

    // Eck et al. (1995). Multiresolution analysis of arbitrary meshes.
    //
    // These constant corresponds to a convex combination of elements with
    // weights equal to the cotangens of the angle opposite to the edge {i, j}.
    // See (Floater, 2003) for the weights formula based on the angle.
    case Default:
    case Harmonic: {
      vtkIdType npts, *pts;
      vtkIdType ptId1 = static_cast<vtkIdType>(i);
      vtkIdType ptId2 = static_cast<vtkIdType>(j);
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
      mirtkAssert(cellIds->GetNumberOfIds() == 2, "surface is triangular mesh");
      double p1[3], p2[3], p3[3], l12, l13, l23;
      _Surface->GetPoint(ptId1, p1);
      _Surface->GetPoint(ptId2, p2);
      l12 = vtkMath::Distance2BetweenPoints(p1, p2);
      for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
        _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
        while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
        _Surface->GetPoint(pts[0], p3);
        l13 = vtkMath::Distance2BetweenPoints(p1, p3);
        l23 = vtkMath::Distance2BetweenPoints(p2, p3);
        k_ij += (l13 + l23 - l12) / vtkTriangle::TriangleArea(p1, p2, p3);
      }
    } break;

    // Wachpress (1975). A Rational Finite Element Basis.
    // Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
    case BarycentricCoordinates: {
      vtkIdType npts, *pts;
      vtkIdType ptId1 = static_cast<vtkIdType>(i);
      vtkIdType ptId2 = static_cast<vtkIdType>(j);
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
      mirtkAssert(cellIds->GetNumberOfIds() == 2, "surface is triangular mesh");
      double p1[3], p2[3], p3[3], dist2;
      _Surface->GetPoint(ptId1, p1);
      _Surface->GetPoint(ptId2, p2);
      dist2 = vtkMath::Distance2BetweenPoints(p1, p2);
      for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
        _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
        while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
        _Surface->GetPoint(pts[0], p3);
        k_ij += Cotangent(p1, p2, p3) / dist2;
      }
    } break;

    // Floater (2003). Mean value coordinates.
    case MeanValueCoordinates: {
      vtkIdType npts, *pts;
      vtkIdType ptId1 = static_cast<vtkIdType>(i);
      vtkIdType ptId2 = static_cast<vtkIdType>(j);
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      _Surface->GetCellEdgeNeighbors(-1, ptId1, ptId2, cellIds);
      mirtkAssert(cellIds->GetNumberOfIds() == 2, "surface is triangular mesh");
      double p1[3], p2[3], p3[3], e1[3], e2[3], l12;
      _Surface->GetPoint(ptId1, p1);
      _Surface->GetPoint(ptId2, p2);
      vtkMath::Subtract(p2, p1, e1);
      l12 = vtkMath::Normalize(e1);
      for (vtkIdType cellIdx = 0; cellIdx < cellIds->GetNumberOfIds(); ++cellIdx) {
        _Surface->GetCellPoints(cellIds->GetId(cellIdx), npts, pts);
        while (pts[0] == ptId1 || pts[0] == ptId2) ++pts;
        _Surface->GetPoint(pts[0], p3);
        vtkMath::Subtract(p3, p1, e2);
        vtkMath::Normalize(e2);
        k_ij += tan(.5 * acos(vtkMath::Dot(e1, e2))) / l12;
      }
    } break;

    // Mentioned in Eck et al. (1995) that these were used by Kent et al.
    case InverseLength: {
      double p1[3], p2[3];
      vtkIdType ptId1 = static_cast<vtkIdType>(i);
      vtkIdType ptId2 = static_cast<vtkIdType>(j);
      _Surface->GetPoint(ptId1, p1);
      _Surface->GetPoint(ptId2, p2);
      k_ij = 1.0 / vtkMath::Distance2BetweenPoints(p1, p2);
    } break;
  }
  return k_ij;
}

// -----------------------------------------------------------------------------
bool LinearSpringSurfaceMapper::Remesh()
{
  if (_SpringConstant == Harmonic || _SpringConstant == MeanValueCoordinates) {
    if (!IsTriangularMesh(_Surface)) {
      _Surface = Triangulate(_Surface);
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
void LinearSpringSurfaceMapper::Solve()
{
  typedef Eigen::MatrixXd                  Values;
  typedef Eigen::SparseMatrix<double>      Matrix;
  typedef Eigen::Triplet<double>           NZEntry;

  const int n = NumberOfFreePoints();
  const int m = NumberOfComponents();

  int i, j, r, c, l;

  const bool symmetric = this->Symmetric();

  Matrix A(n, n);
  Values b(n, m);
  {
    EdgeTable edgeTable(_Surface);
    EdgeIterator edgeIt(edgeTable);

    b.setZero();
    Array<NZEntry> k;
    Array<double>  k_ii(n, .0);
    k.reserve(2 * edgeTable.NumberOfEdges() + n);

    double k_ij, k_ji;
    for (edgeIt.InitTraversal(); edgeIt.GetNextEdge(i, j) != -1;) {
      r = FreePointIndex(i);
      c = FreePointIndex(j);
      if (r >= 0 || c >= 0) {
        k_ij = K(i, j);
        if (symmetric) k_ji = k_ij;
        else           k_ji = K(j, i);
        if (r >= 0 && c >= 0) {
          k.push_back(NZEntry(r, c, -k_ij));
          k.push_back(NZEntry(c, r, -k_ji));
        } else if (r >= 0) {
          for (l = 0; l < m; ++l) {
            b(r, l) += k_ij * GetValue(j, l);
          }
        } else if (c >= 0) {
          for (l = 0; l < m; ++l) {
            b(c, l) += k_ji * GetValue(i, l);
          }
        }
        if (r >= 0) {
          k_ii[r] += k_ij;
        }
        if (c >= 0) {
          k_ii[c] += k_ji;
        }
      }
    }
    for (r = 0; r < n; ++r) {
      k.push_back(NZEntry(r, r, k_ii[r]));
    }

    A.setFromTriplets(k.begin(), k.end());
  }

  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points             = " << _Surface->GetNumberOfPoints() << "\n";
    cout << "  No. of free points                = " << n << "\n";
    cout << "  No. of non-zero stiffness values  = " << A.nonZeros() << "\n";
    cout << "  Dimension of surface map codomain = " << m << "\n";
    cout.flush();
  }

  Values x(n, m);
  for (r = 0; r < n; ++r) {
    i = FreePointId(r);
    for (l = 0; l < m; ++l) {
      x(r, l) = GetValue(i, l);
    }
  }

  int    niter = 0;
  double error = .0;

  if (symmetric) {
    Eigen::ConjugateGradient<Matrix> solver(A);
    if (_NumberOfIterations >  0) solver.setMaxIterations(_NumberOfIterations);
    if (_Tolerance          > .0) solver.setTolerance(_Tolerance);
    x = solver.solveWithGuess(b, x);
    niter = solver.iterations();
    error = solver.error();
  } else {
    Eigen::BiCGSTAB<Matrix> solver(A);
    if (_NumberOfIterations >  0) solver.setMaxIterations(_NumberOfIterations);
    if (_Tolerance          > .0) solver.setTolerance(_Tolerance);
    x = solver.solveWithGuess(b, x);
    niter = solver.iterations();
    error = solver.error();
  }

  if (verbose) {
    cout << "  No. of iterations                 = " << niter << "\n";
    cout << "  Estimated error                   = " << error << "\n";
    cout.flush();
  }

  for (r = 0; r < n; ++r) {
    i = FreePointId(r);
    for (l = 0; l < m; ++l) {
      SetValue(i, l, x(r, l));
    }
  }
}


} // namespace mirtk

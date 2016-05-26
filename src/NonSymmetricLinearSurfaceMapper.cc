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

#include "mirtk/NonSymmetricLinearSurfaceMapper.h"

#include "mirtk/EdgeTable.h"

#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void NonSymmetricLinearSurfaceMapper::CopyAttributes(const NonSymmetricLinearSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
NonSymmetricLinearSurfaceMapper::NonSymmetricLinearSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
NonSymmetricLinearSurfaceMapper::NonSymmetricLinearSurfaceMapper(const NonSymmetricLinearSurfaceMapper &other)
:
  LinearSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NonSymmetricLinearSurfaceMapper &NonSymmetricLinearSurfaceMapper::operator =(const NonSymmetricLinearSurfaceMapper &other)
{
  if (this != &other) {
    LinearSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NonSymmetricLinearSurfaceMapper::~NonSymmetricLinearSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void NonSymmetricLinearSurfaceMapper::Weights(int i, const int *j, double *w_i, int d_i) const
{
  for (int k = 0; k < d_i; ++k) {
    w_i[k] = this->Weight(i, j[k]);
  }
}

// -----------------------------------------------------------------------------
double NonSymmetricLinearSurfaceMapper::Weight(int i, int j) const
{
  cerr << this->NameOfClass() << "::Weight[s]: Either one of the overloads must be implemented in subclass" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void NonSymmetricLinearSurfaceMapper::Solve()
{
  typedef Eigen::MatrixXd             Values;
  typedef Eigen::SparseMatrix<double> Matrix;
  typedef Eigen::Triplet<double>      NZEntry;

  const int n = NumberOfFreePoints();
  const int m = NumberOfComponents();

  int i, k, l, r, c;
  const int *j;

  Matrix A(n, n);
  Values b(n, m);
  {
    EdgeTable edgeTable(_Surface);
    EdgeIterator edgeIt(edgeTable);

    b.setZero();
    Array<NZEntry> w;
    Array<double>  w_ii(n, .0);
    w.reserve(2 * edgeTable.NumberOfEdges() + n);

    int     d_i;
    double *w_i = new double[edgeTable.MaxNumberOfAdjacentPoints()];
    for (i = 0; i < NumberOfPoints(); ++i) {
      r = FreePointIndex(i);
      if (r >= 0) {
        edgeTable.GetAdjacentPoints(i, d_i, j);
        this->Weights(i, j, w_i, d_i);
        for (k = 0; k < d_i; ++k) {
          c = FreePointIndex(j[k]);
          if (c >= 0) {
            w.push_back(NZEntry(r, c, -w_i[k]));
          } else {
            for (l = 0; l < m; ++l) {
              b(r, l) += w_i[k] * GetValue(j[k], l);
            }
          }
          w_ii[r] += w_i[k];
        }
      }
    }
    for (r = 0; r < n; ++r) {
      w.push_back(NZEntry(r, r, w_ii[r]));
    }
    delete[] w_i;

    A.setFromTriplets(w.begin(), w.end());
  }

  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points             = " << NumberOfPoints() << "\n";
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

  Eigen::BiCGSTAB<Matrix> solver(A);
  if (_NumberOfIterations >  0) solver.setMaxIterations(_NumberOfIterations);
  if (_Tolerance          > .0) solver.setTolerance(_Tolerance);
  x = solver.solveWithGuess(b, x);
  niter = solver.iterations();
  error = solver.error();

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

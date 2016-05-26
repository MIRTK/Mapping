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

#include "mirtk/ConformalFlatteningSurfaceMapper.h"

#include "mirtk/VtkMath.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PolyDataCurvature.h"

#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper
::CopyAttributes(const ConformalFlatteningSurfaceMapper &other)
{
  _PolarCellId = other._PolarCellId;
  _MapToSphere = other._MapToSphere;
  _Scale       = other._Scale;
  _Radius      = other._Radius;
}

// -----------------------------------------------------------------------------
ConformalFlatteningSurfaceMapper::ConformalFlatteningSurfaceMapper()
:
  _PolarCellId(-1),
  _MapToSphere(true),
  _Scale(.0),
  _Radius(1.0)
{
}

// -----------------------------------------------------------------------------
ConformalFlatteningSurfaceMapper::ConformalFlatteningSurfaceMapper(
  const ConformalFlatteningSurfaceMapper &other
) :
  LinearSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ConformalFlatteningSurfaceMapper &ConformalFlatteningSurfaceMapper
::operator =(const ConformalFlatteningSurfaceMapper &other)
{
  if (this != &other) {
    LinearSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ConformalFlatteningSurfaceMapper::~ConformalFlatteningSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper::Initialize()
{
  // Initialize base class
  LinearSurfaceMapper::Initialize();

  // Check that fixed point cell ID is valid
  if (_PolarCellId >= _Surface->GetNumberOfCells()) {
    cerr << this->NameOfType() << "::Initialize: Invalid fixed point cell ID!" << endl;
    exit(1);
  }

  // Input surface must have genus 0
  if (Genus(_Surface) != .0) {
    cerr << this->NameOfType() << "::Initialize: Input surface must have genus 0!" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper::InitializeValues()
{
  if (_Input) {
    _Values = _Input->NewInstance();
    _Values->SetName(_Input->GetName());
  } else {
    _Values = vtkSmartPointer<vtkFloatArray>::New();
    _Values->SetName("Map");
  }
  _Values->SetNumberOfComponents(_MapToSphere ? 3 : 2);
  _Values->SetNumberOfTuples(_Surface->GetNumberOfPoints());
  for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
    _Values->FillComponent(j, .0);
  }
}

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper::InitializeMask()
{
  // Choose cell containing polar point automatically
  if (_PolarCellId < 0) {
    PolyDataCurvature filter(PolyDataCurvature::Mean);
    filter.Input(_Surface);
    filter.Run();
    vtkDataArray *curvature = filter.GetMeanCurvature();
    for (vtkIdType ptId = 1; ptId < curvature->GetNumberOfTuples(); ++ptId) {
      curvature->SetComponent(ptId, 0, abs(curvature->GetComponent(ptId, 0)));
    }
    vtkIdType polarPoint = 0;
    double    polarValue = curvature->GetComponent(0, 0);
    for (vtkIdType ptId = 1; ptId < curvature->GetNumberOfTuples(); ++ptId) {
      if (curvature->GetComponent(ptId, 0) < polarValue) {
        polarValue = curvature->GetComponent(ptId, 0);
        polarPoint = ptId;
      }
    }
    unsigned short ncells;
    vtkIdType      npts, *cells, *pts;
    _Surface->GetPointCells(polarPoint, ncells, cells);
    _PolarCellId = cells[0];
    double avgValue, minValue = .0;
    for (unsigned short i = 0; i < ncells; ++i) {
      _Surface->GetCellPoints(cells[i], npts, pts);
      avgValue = .0;
      for (vtkIdType j = 0; j < npts; ++j) {
        avgValue += curvature->GetComponent(pts[j], 0);
      }
      avgValue /= npts;
      if (i == 0 || avgValue < minValue) {
        _PolarCellId = cells[i];
        minValue     = avgValue;
      }
    }
  }

  // Mark corners of cell with polar point as points with fixed values
  vtkIdType npts, *pts;
  _Surface->GetCellPoints(_PolarCellId, npts, pts);

  _Fixed = vtkSmartPointer<vtkUnsignedCharArray>::New();
  _Fixed->SetName("FixedPoints");
  _Fixed->SetNumberOfComponents(1);
  _Fixed->SetNumberOfTuples(_Surface->GetNumberOfPoints());

  _Fixed->FillComponent(0, .0);
  for (vtkIdType i = 0; i < npts; ++i) {
    _Fixed->SetComponent(pts[i], 0, 1.0);
  }
}

// -----------------------------------------------------------------------------
bool ConformalFlatteningSurfaceMapper::Remesh()
{
  // Triangulate surface if necessary
  if (!IsTriangularMesh(_Surface)) {
    _Surface = Triangulate(_Surface);
    if (_PolarCellId >= 0) {
      // TODO: Map _PolarCellId to corresponding triangle ID
      cerr << this->NameOfType() << "::Remesh: Polar cell ID to triangle ID of new mesh conversion not implemented!" << endl;
      cerr << "  Triangulate surface mesh before running this filter." << endl;
      exit(1);
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper::Solve()
{
  typedef Eigen::MatrixXd             Values;
  typedef Eigen::SparseMatrix<double> Matrix;

  const int n = NumberOfPoints();
  const int m = 2;

  // Calculate matrix D
  Matrix D(n, n);
  {
    typedef Eigen::VectorXi Vector;

    vtkIdType npts, *pts;                // cell links list reference
    double posA [3], posB [3], posC [3]; // position of cell corner points
    double ctgABC, ctgBCA, ctgCAB;       // cotangent of cell angles

    {
      EdgeTable edgeTable(_Surface);
      D.reserve(Vector::Constant(n, edgeTable.MaxNumberOfAdjacentPoints() + 1));
    }
    for (vtkIdType cellId = 0; cellId < _Surface->GetNumberOfCells(); ++cellId) {
      _Surface->GetCellPoints(cellId, npts, pts);
      const vtkIdType &ptIdA = pts[0];
      const vtkIdType &ptIdB = pts[1];
      const vtkIdType &ptIdC = pts[2];

      _Surface->GetPoint(ptIdA, posA);
      _Surface->GetPoint(ptIdB, posB);
      _Surface->GetPoint(ptIdC, posC);

      ctgABC = Cotangent(posA, posB, posC);
      ctgBCA = Cotangent(posB, posC, posA);
      ctgCAB = Cotangent(posC, posA, posB);

      D.coeffRef(ptIdA, ptIdA) += ctgABC + ctgBCA;
      D.coeffRef(ptIdA, ptIdB) -= ctgBCA;
      D.coeffRef(ptIdA, ptIdC) -= ctgABC;

      D.coeffRef(ptIdB, ptIdB) += ctgBCA + ctgCAB;
      D.coeffRef(ptIdB, ptIdA) -= ctgBCA;
      D.coeffRef(ptIdB, ptIdC) -= ctgCAB;

      D.coeffRef(ptIdC, ptIdC) += ctgCAB + ctgABC;
      D.coeffRef(ptIdC, ptIdB) -= ctgCAB;
      D.coeffRef(ptIdC, ptIdA) -= ctgABC;
    }
    D.makeCompressed();
  }

  // Calculate (complex) right hand side vector b
  Values b(n, m);
  {
    vtkIdType npts, *pts;             // polar cell links list reference
    double posA[3], posB[3], posC[3]; // position of cell corner points
    double vecAB[3], vecCA[3];        // edge vectors AB and AC
    double posE[3];                   // projection of C onto AB
    double vecEC[3];                  // vector from E to C

    // Get corners of cell containing polar point
    _Surface->GetCellPoints(_PolarCellId, npts, pts);
    const vtkIdType &ptIdA = pts[0];
    const vtkIdType &ptIdB = pts[1];
    const vtkIdType &ptIdC = pts[2];

    _Surface->GetPoint(ptIdA, posA);
    _Surface->GetPoint(ptIdB, posB);
    _Surface->GetPoint(ptIdC, posC);

    // Edge vectors AB and AC
    vtkMath::Subtract(posB, posA, vecAB);
    vtkMath::Subtract(posA, posC, vecCA);

    // Orthogonal projection of C onto AB
    const double normAB2 = vtkMath::Dot(vecAB, vecAB);
    const double theta   = - vtkMath::Dot(vecCA, vecAB) / normAB2;
    posE[0] = posA[0] + theta * vecAB[0];
    posE[1] = posA[1] + theta * vecAB[1];
    posE[2] = posA[2] + theta * vecAB[2];
    vtkMath::Subtract(posC, posE, vecEC);

    // Constraints of corner points
    b.setZero();

    const double x = 2.0 / sqrt(normAB2);
    b(ptIdA, 0) = -x;
    b(ptIdB, 0) =  x;

    const double y = 2.0 / vtkMath::Norm(vecEC);
    b(ptIdA, 1) =  y * (1.0 - theta);
    b(ptIdB, 1) =  y * theta;
    b(ptIdC, 1) = -y;
  }

  // Solve linear system
  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points  = " << NumberOfPoints() << "\n";
    cout << "  No. of fixed points    = " << NumberOfFixedPoints() << "\n";
    cout << "  No. of non-zero values = " << D.nonZeros() << "\n";
    cout << "  Dimension of codomain  = " << m << "\n";
    cout.flush();
  }

  int    niter = 0;
  double error = .0;

  Eigen::ConjugateGradient<Matrix> solver(D);
  if (_NumberOfIterations >  0) solver.setMaxIterations(_NumberOfIterations);
  if (_Tolerance          > .0) solver.setTolerance(_Tolerance);
  Values x = solver.solve(b);
  niter = solver.iterations();
  error = solver.error();

  if (verbose) {
    cout << "  No. of iterations      = " << niter << "\n";
    cout << "  Estimated error        = " << error << "\n";
    cout.flush();
  }

  for (int i = 0; i < n; ++i) {
    for (int l = 0; l < m; ++l) {
      SetValue(i, l, x(i, l));
    }
  }
}

// -----------------------------------------------------------------------------
void ConformalFlatteningSurfaceMapper::Finalize()
{
  if (_MapToSphere) {
    const int n = NumberOfPoints();
    double x, y, r2, scale = _Scale;

    // Choose scale for inverse stereographic projection automatically s.t.
    // afterwards, upper and lower hemisphere have same number of points
    if (scale <= .0) {
      Array<double> v_r2(n);
      auto v_r2_it = v_r2.begin();
      for (int i = 0; i < n; ++i, ++v_r2_it) {
        x = GetValue(i, 0);
        y = GetValue(i, 1);
        *v_r2_it = x*x + y*y;
      }
      sort(v_r2.begin(), v_r2.end());
      int i = ((n % 2) == 0 ? n : n - 1) / 2;
      scale = 1.0 / sqrt(v_r2[i]);
    }

    // Perform inverse stereographic projection
    const double s = 2.0 * _Radius;
    for (int i = 0; i < n; ++i) {
      x = scale * GetValue(i, 0);
      y = scale * GetValue(i, 1);
      r2 = x*x + y*y;
      SetValue(i, 0, s * x / (1.0 + r2));
      SetValue(i, 1, s * y / (1.0 + r2));
      SetValue(i, 2, s * r2 / ( 1.0 + r2) - 1.0);
    }
  }

  // Finalize base class
  LinearSurfaceMapper::Finalize();
}


} // namespace mirtk

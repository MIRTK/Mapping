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

#include "mirtk/NearOptimalSurfaceMapper.h"

#include "mirtk/IntrinsicSurfaceMapper.h"

#include "vtkPointData.h"
#include "vtkCellData.h"

#if MIRTK_USE_FLOAT_BY_DEFAULT
  #include "vtkFloatArray.h"
#else
  #include "vtkDoubleArray.h"
#endif


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void NearOptimalSurfaceMapper::CopyAttributes(const NearOptimalSurfaceMapper &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;

  if (other._Values) {
    _Values = other._Values->NewInstance();
    _Values->DeepCopy(other._Values);
  } else {
    _Values = nullptr;
  }
}

// -----------------------------------------------------------------------------
NearOptimalSurfaceMapper::NearOptimalSurfaceMapper()
:
  _NumberOfIterations(0),
  _Tolerance(-1.)
{
}

// -----------------------------------------------------------------------------
NearOptimalSurfaceMapper::NearOptimalSurfaceMapper(const NearOptimalSurfaceMapper &other)
:
  FixedBoundarySurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NearOptimalSurfaceMapper &NearOptimalSurfaceMapper
::operator =(const NearOptimalSurfaceMapper &other)
{
  if (this != &other) {
    FixedBoundarySurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NearOptimalSurfaceMapper::~NearOptimalSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void NearOptimalSurfaceMapper::Initialize()
{
  // Initialize base class
  FixedBoundarySurfaceMapper::Initialize();

  // Check input
  if (_Surface->GetNumberOfPolys() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input point set must be a surface mesh" << endl;
    exit(1);
  }

  // Initialize map values and
  // determine sets of points with free or fixed map values, respectively
  const int m = NumberOfComponents();
  const int n = NumberOfPoints();

  #if MIRTK_USE_FLOAT_BY_DEFAULT
    _Values = vtkSmartPointer<vtkFloatArray>::New();
  #else
    _Values = vtkSmartPointer<vtkDoubleArray>::New();
  #endif
  _Values->SetName("SurfaceMap");
  _Values->SetNumberOfComponents(m);
  _Values->SetNumberOfTuples(static_cast<vtkIdType>(n));
}

// -----------------------------------------------------------------------------
void NearOptimalSurfaceMapper::ComputeMap()
{
  const int       m = NumberOfComponents();
  const vtkIdType n = NumberOfPoints();

  vtkSmartPointer<vtkDataArray> u = _Values;
  vtkSmartPointer<vtkDataArray> v = _Values->NewInstance();
  v->SetNumberOfComponents(m);
  v->SetNumberOfTuples(n);

  // Compute intrinsic surface maps
  {
    IntrinsicSurfaceMapper mapper;
    mapper.NumberOfIterations(_NumberOfIterations);
    mapper.Tolerance(_Tolerance);
    mapper.Surface(_Surface);
    mapper.EdgeTable(_EdgeTable);
    mapper.Boundary(_Boundary);
    mapper.Input(_Input);
    mapper.Initialize();

    // Compute discrete conformal map
    mapper.Lambda(1.);
    mapper.ComputeMap();
    for (vtkIdType i = 0; i < n; ++i)
    for (int       j = 0; j < m; ++j) {
      u->SetComponent(i, j, mapper.GetValue(static_cast<int>(i), j));
    }

    // Compute discrete authalic map
    mapper.Lambda(0.);
    mapper.ComputeMap();
    for (vtkIdType i = 0; i < n; ++i)
    for (int       j = 0; j < m; ++j) {
      v->SetComponent(i, j, mapper.GetValue(static_cast<int>(i), j));
    }
  }

  // Determine optimal lambda
  const double lambda = this->ComputeLambda(u, v);
  const double mu     = 1. - lambda;

  // Set near-optimal surface map values
  double value;
  for (vtkIdType i = 0; i < n; ++i)
  for (int       j = 0; j < m; ++j) {
    value = lambda * u->GetComponent(i, j) + mu * v->GetComponent(i, j);
    _Values->SetComponent(i, j, value);
  }
}

// -----------------------------------------------------------------------------
void NearOptimalSurfaceMapper::Finalize()
{
  SharedPtr<PiecewiseLinearMap> map = NewShared<PiecewiseLinearMap>();
  vtkSmartPointer<vtkPolyData> domain = _Surface->NewInstance();
  domain->ShallowCopy(_Surface);
  domain->GetPointData()->Initialize();
  domain->GetCellData()->Initialize();
  map->Domain(domain);
  map->Values(_Values);
  _Output = map;
}


} // namespace mirtk

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

#include "mirtk/UniformBoundaryParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
UniformBoundaryParameterizer::UniformBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
UniformBoundaryParameterizer
::UniformBoundaryParameterizer(const UniformBoundaryParameterizer &other)
:
  BoundaryParameterizer(other)
{
}

// -----------------------------------------------------------------------------
UniformBoundaryParameterizer &UniformBoundaryParameterizer
::operator =(const UniformBoundaryParameterizer &other)
{
  if (this != &other) {
    BoundaryParameterizer::operator =(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
UniformBoundaryParameterizer::~UniformBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundaryParameterizer *UniformBoundaryParameterizer::NewCopy() const
{
  return new UniformBoundaryParameterizer(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void UniformBoundaryParameterizer::Parameterize()
{
  const int npoints   = NumberOfBoundaryPoints();
  const int nselected = NumberOfSelectedPoints();
  const int i0        = (nselected > 0 ? SelectedPointIndex(0) : 0);

  // Uniform parameterization of piecewise linear curve
  double t = .0, dt = 1.0 / npoints;
  for (int n = 0, i = i0; n < npoints; ++n, t += dt) {
    _Values[i] = t;
    if (++i == npoints) i = 0;
  }

  // Revert orientation of curve if points where selected in reverse order
  if (nselected > 2) {
    double t1 = _Values[SelectedPointIndex(1)];
    double t2 = _Values[SelectedPointIndex(2)];
    if (t2 < t1) {
      for (int n = 0; n < npoints; ++n) {
        if (n != i0) _Values[n] = 1.0 - _Values[n];
      }
    }
  }
}


} // namespace mirtk

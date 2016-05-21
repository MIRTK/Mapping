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

#include "mirtk/LinearBoundaryParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LinearBoundaryParameterizer::LinearBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
LinearBoundaryParameterizer
::LinearBoundaryParameterizer(const LinearBoundaryParameterizer &other)
:
  BoundaryParameterizer(other)
{
}

// -----------------------------------------------------------------------------
LinearBoundaryParameterizer &LinearBoundaryParameterizer
::operator =(const LinearBoundaryParameterizer &other)
{
  if (this != &other) {
    BoundaryParameterizer::operator =(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LinearBoundaryParameterizer::~LinearBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundaryParameterizer *LinearBoundaryParameterizer::NewCopy() const
{
  return new LinearBoundaryParameterizer(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void LinearBoundaryParameterizer::Parameterize()
{
  const int npoints   = NumberOfBoundaryPoints();
  const int nselected = NumberOfSelectedPoints();
  const int nfixed    = (nselected > 0 ? nselected             : 1);
  const int i0        = (nselected > 0 ? SelectedPointIndex(0) : 0);

  // Compute edge lengths
  const Vector l = BoundaryEdgeLengths();

  // Compute distance of each point from point with index i0
  double d = .0;
  Vector d1(nfixed);
  for (int n = 0, s = 0, i = i0; n < npoints; ++n) {
    _Values[i] = d, d += l(i);
    if (++i == npoints) i = 0;
    if (i == i0 || IsSelected(i)) {
      d1(s) = d, ++s;
    }
  }

  // Map distances to curve parameter values using the linear function:
  //
  //   t(d) = ((d - d0) / (d1 - d0)) * (t1 - t0) + t0
  //
  // where d0 and d1 are the distances of the nearest selected points,
  // and t0 and t1 are the respective parameter values of these points.
  double tdiff = 1.0 / nfixed;
  double d0    = .0;
  double scale = tdiff / d1(0);
  for (int n = 0, s = 0, i = i0; n < npoints; ++n) {
    _Values[i] = scale * (_Values[i] - d0) + s * tdiff;
    if (++i == npoints) i = 0;
    if (i != i0 && IsSelected(i)) {
      d0 = d1(s), ++s;
      scale = tdiff / (d1(s) - d0);
    }
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

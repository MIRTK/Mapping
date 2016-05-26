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

#include "mirtk/ChordLengthBoundaryParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ChordLengthBoundaryParameterizer::ChordLengthBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
ChordLengthBoundaryParameterizer
::ChordLengthBoundaryParameterizer(const ChordLengthBoundaryParameterizer &other)
:
  BoundaryParameterizer(other)
{
}

// -----------------------------------------------------------------------------
ChordLengthBoundaryParameterizer &ChordLengthBoundaryParameterizer
::operator =(const ChordLengthBoundaryParameterizer &other)
{
  if (this != &other) {
    BoundaryParameterizer::operator =(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ChordLengthBoundaryParameterizer::~ChordLengthBoundaryParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundaryParameterizer *ChordLengthBoundaryParameterizer::NewCopy() const
{
  return new ChordLengthBoundaryParameterizer(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ChordLengthBoundaryParameterizer::Parameterize()
{
  const int npoints   = NumberOfBoundaryPoints();
  const int nselected = NumberOfSelectedPoints();

  // Compute edge lengths
  Vector l = BoundaryEdgeLengths();

  // Chord length parameterization
  double t = .0;
  const double L = l.Sum();
  for (int i = 0; i < npoints; ++i) {
    _Values[i] = t;
    t += l(i) / L;
  }

  // Reparameterize such that first selected point has t=0
  if (nselected > 1) {
    const int    i0 = SelectedPointIndex(0);
    const double t0 = _Values[i0];
    if (t0 != .0) {
      for (int i = 0; i < npoints; ++i) {
        _Values[i] -= t0;
        if (_Values[i] < .0) _Values[i] += 1.0;
      }
      _Values[i0] = .0;
    }
  }

  // Revert orientation of curve if points where selected in reverse order
  if (nselected > 2) {
    double t1 = _Values[SelectedPointIndex(1)];
    double t2 = _Values[SelectedPointIndex(2)];
    if (t2 < t1) {
      for (int n = 0; n < npoints; ++n) {
        if (_Values[n] != .0) {
          _Values[n] = 1.0 - _Values[n];
        }
      }
    }
  }
}


} // namespace mirtk

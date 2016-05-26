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

#include "mirtk/LinearSurfaceMapper.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void LinearSurfaceMapper::CopyAttributes(const LinearSurfaceMapper &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
}

// -----------------------------------------------------------------------------
LinearSurfaceMapper::LinearSurfaceMapper()
:
  _NumberOfIterations(0),
  _Tolerance(.0)
{
}

// -----------------------------------------------------------------------------
LinearSurfaceMapper::LinearSurfaceMapper(const LinearSurfaceMapper &other)
:
  SurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LinearSurfaceMapper &LinearSurfaceMapper::operator =(const LinearSurfaceMapper &other)
{
  if (this != &other) {
    SurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LinearSurfaceMapper::~LinearSurfaceMapper()
{
}


} // namespace mirtk

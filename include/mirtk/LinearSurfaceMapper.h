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

#ifndef MIRTK_LinearSurfaceMapper_H
#define MIRTK_LinearSurfaceMapper_H

#include "mirtk/SurfaceMapper.h"


namespace mirtk {


/**
 * Obtains a surface map as the solution of a sparse linear system of equations
 *
 * For a given choice of edge weights, a sparse system of linear equations is solved
 * iteratively to obtain a surface parameterization (Marchandise et al., 2014).
 */
class LinearSurfaceMapper : public SurfaceMapper
{
  mirtkAbstractMacro(LinearSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Tolerance for sparse linear solver
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LinearSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  LinearSurfaceMapper();

  /// Copy constructor
  LinearSurfaceMapper(const LinearSurfaceMapper &);

  /// Assignment operator
  LinearSurfaceMapper &operator =(const LinearSurfaceMapper &);

public:

  /// Destructor
  virtual ~LinearSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Construct and solve linear system of equations
  virtual void Solve() = 0;

};


} // namespace mirtk

#endif // MIRTK_LinearSurfaceMapper_H

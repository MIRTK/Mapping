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

#ifndef MIRTK_LeastSquaresConformalSurfaceMapper_H
#define MIRTK_LeastSquaresConformalSurfaceMapper_H

#include "mirtk/SymmetricLinearSurfaceMapper.h"


namespace mirtk {


/**
 * Levy's least squares conformal mapping method
 *
 * The mapping method implemented by this class (Levy, 2002) is related to the
 * harmonic mapping methods by Pinkall and Polthier (1993), Eck et al. (1995),
 * and Haker et al. (2000). Unlike these authors, Levy et al. consider the
 * conformality condition on each triangle and link the (u, v) parameters of
 * the resulting map directly rather than indirectly via the right-hand side
 * of the systems of linear equations resulting from the discretization of the
 * Laplace-Beltrami operator at the vertices using the finite element method (FEM).
 * Moreover, Levy et al. add the constraints that edges should map to straight
 * lines in which case conformality cannot be guaranteed. Their method finds
 * instead an as-conformal-as-possible solution in the least squares sense.
 *
 * - Lévy et al. (2002). Least squares conformal maps for automatic texture atlas
 *   generation. ACM Trans. Graphics, 21(3), 362–371.
 * - Pinkall and Polthier (1993). See IntrinsicParameterizationSurfaceMapper.
 * - Eck et al. (1995). See IntrinsicParameterizationSurfaceMapper.
 * - Haker et al. (2000). See ConformalFlatteningSurfaceMapper class.
 *
 * \sa IntrinsicParameterizationSurfaceMapper, ConformalFlatteningSurfaceMapper
 */
class LeastSquaresConformalSurfaceMapper : public SymmetricLinearSurfaceMapper
{
  mirtkObjectMacro(LeastSquaresConformalSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LeastSquaresConformalSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  LeastSquaresConformalSurfaceMapper();

  /// Copy constructor
  LeastSquaresConformalSurfaceMapper(const LeastSquaresConformalSurfaceMapper &);

  /// Assignment operator
  LeastSquaresConformalSurfaceMapper &operator =(const LeastSquaresConformalSurfaceMapper &);

  /// Destructor
  virtual ~LeastSquaresConformalSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of undirected edge (i, j)
  ///
  /// \param[in] i First end point.
  /// \param[in] j Second end point.
  ///
  /// \returns Weight of undirected edge (i, j).
  virtual double Weight(int i, int j) const;

};


} // namespace mirtk

#endif // MIRTK_LeastSquaresConformalSurfaceMapper_H

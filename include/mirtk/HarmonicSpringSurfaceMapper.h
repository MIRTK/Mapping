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

#ifndef MIRTK_HarmonicSpringSurfaceMapper_H
#define MIRTK_HarmonicSpringSurfaceMapper_H

#include "mirtk/SymmetricLinearSurfaceMapper.h"


namespace mirtk {


/**
 * Harmonic surface map as the critical point of a linear spring network energy
 *
 * This surface mapper finds a piecewise linear approximation of the harmonic
 * map of a non-closed surface given a fixed boundary map, i.e., with Dirichlet
 * boundary conditions. The linear spring model energy function corresponds
 * to a discretization of the Dirichlet energy functional and the weights
 * are identical to those obtained by the classic finite element method (FEM).
 * The cotangent weights are numerically differently computed by this mapper
 * according to the spring constants presented in Eck et al. (1995). In general,
 * the IntrinsicSurfaceParameterization implementation should be preferred
 * as it is more generic.
 *
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 *
 * @sa IntrinsicSurfaceParameterization
 */
class HarmonicSpringSurfaceMapper : public SymmetricLinearSurfaceMapper
{
  mirtkObjectMacro(HarmonicSpringSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const HarmonicSpringSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  HarmonicSpringSurfaceMapper();

  /// Copy constructor
  HarmonicSpringSurfaceMapper(const HarmonicSpringSurfaceMapper &);

  /// Assignment operator
  HarmonicSpringSurfaceMapper &operator =(const HarmonicSpringSurfaceMapper &);

public:

  /// Destructor
  virtual ~HarmonicSpringSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Remesh surface if necessary
  virtual bool Remesh();

  /// Weight of undirected edge (i, j)
  ///
  /// \param[in] i First end point.
  /// \param[in] j Second end point.
  ///
  /// \returns Weight of undirected edge (i, j).
  virtual double Weight(int i, int j) const;

};


} // namespace mirtk

#endif // MIRTK_HarmonicSpringSurfaceMapper_H

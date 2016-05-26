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

#ifndef MIRTK_ConformalFlatteningSurfaceMapper_H
#define MIRTK_ConformalFlatteningSurfaceMapper_H

#include "mirtk/LinearSurfaceMapper.h"


namespace mirtk {


/**
 * Computes discrete conformal map of (closed) surface to the complex plane
 *
 * The conformal map is the solution of a sparse system of linear equations,
 * where the edge weights are identical to those used by Pinkall and Polthier
 * (1993) and later Eck et al. (1995). These cotangent weights are obtained
 * using the finite element method (FEM) to discretize the Laplace-Beltrami
 * operator. The main difference between the harmonic mapping method of Pinkall
 * and Polthier and the one implemented by this mapper (Angenent & Haker, 1999-2000)
 * are the boundary constraints. While Pinkall and Polthier, as well as Eck et al.,
 * used Dirichlet boundary conditions on the boundary of the non-closed surface,
 * Angenent and Haker derive the linear system from special boundary constraints
 * intended to map a closed surface to the complex plane and subsequently to the
 * sphere. Here, a single point on the surface is constraint to map to the
 * (north/south) pole of the sphere, i.e., the point for which the stereographic
 * projection is undefined. Dirichlet boundary conditions are applied only to the
 * three vertices of the triangle containing the selected point.
 *
 * The resulting map to the complex plane can be composed with an inverse
 * stereographic projection to obtain a conformal map to the sphere.
 *
 * - Angenent et al. (1999). Conformal geometry and brain flattening. MICCAI, 271–278.
 * - Angenent et al. (n.d.). On the Laplace-Beltrami Operator and Brain Surface Flattening.
 * - Angenent et al. (2000). On Area Preserving Mappings of Minimal Distortion.
 *   In System Theory: Modeling, Analysis and Control, pp. 275–286.
 * - Haker et al. (2000). Conformal surface parameterization for texture mapping.
 *   IEEE Trans. Vis. Comput. Graphics, 6(2), 181–189.
 * - Pinkall and Polthier (1993). Computing Discrete Minimal Surfaces and Their Conjugates.
 *   Experiment. Math., 2(1), 15–36.
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 */
class ConformalFlatteningSurfaceMapper : public LinearSurfaceMapper
{
  mirtkObjectMacro(ConformalFlatteningSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// ID of cell containing polar point P used for inverse stereographic projection
  mirtkPublicAttributeMacro(int, PolarCellId);

  /// Whether to compose flattening map with inverse stereographic projection
  mirtkPublicAttributeMacro(bool, MapToSphere);

  /// Scaling factor applied to complex plane coordinates before projection to sphere
  mirtkPublicAttributeMacro(double, Scale);

  /// Radius of sphere to which flattened map is projected
  mirtkPublicAttributeMacro(double, Radius);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ConformalFlatteningSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ConformalFlatteningSurfaceMapper();

  /// Copy constructor
  ConformalFlatteningSurfaceMapper(const ConformalFlatteningSurfaceMapper &);

  /// Assignment operator
  ConformalFlatteningSurfaceMapper &operator =(const ConformalFlatteningSurfaceMapper &);

  /// Destructor
  virtual ~ConformalFlatteningSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Initialize map values at surface points
  virtual void InitializeValues();

  /// Initialize mask of surface points with fixed map values
  virtual void InitializeMask();

  /// Remesh surface if necessary
  virtual bool Remesh();

  /// Construct and solve symmetric system of linear equations
  virtual void Solve();

  /// Finalize filter execution
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_ConformalFlatteningSurfaceMapper_H

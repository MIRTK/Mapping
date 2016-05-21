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

#ifndef MIRTK_LinearSpringSurfaceMapper_H
#define MIRTK_LinearSpringSurfaceMapper_H

#include "mirtk/SurfaceMapper.h"


namespace mirtk {


/**
 * Computes a surface map using a linear spring model
 *
 * This surface mapper computes a surface map using a linear spring model of the
 * surface mesh of a given piecewise linear complex (PLC). It implements different
 * spring constants proposed in the literature including spring constants derived
 * from the application of the finite element method (FEM) to find an approximate
 * solution to the Laplace equation (a minimizer of the Dirichlet energy functional;
 * Eck et al. 1995), and spring constants which are derived from the convex
 * combination weights resulting from the generalization of barycentric coordinates
 * used to interpolate values within a k-sided polygon (Floater, 2003).
 *
 * For a given choice of spring constants, a sparse system of linear equations is
 * solved iteratively to obtain a surface parameterization (Marchandise et al., 2014).
 *
 * References:
 * - Marchandise et al. (2014). Optimal parametrizations for surface remeshing.
 *   Engineering with Computers, 30(3), 383–402.
 * - Floater (2003). Mean value coordinates. Computer Aided Geometric Design, 20(1):19–37.
 * - Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
 *   Journal of Graphics Tools, 7(1), 13–22.
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 * - Wachpress (1975). A Rational Finite Element Basis.
 * - Tutte (1963). How to draw a graph. Proc. London Math. Soc, 8(May 1962), 743–767.
 */
class LinearSpringSurfaceMapper : public SurfaceMapper
{
  mirtkObjectMacro(LinearSpringSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Enumeration of available spring constants
  enum SpringConstantType
  {
    Default,                ///< Use default spring constants
    Uniform,                ///< Uniform spring constants. The resulting
                            ///< surface map is known as Tutte's embedding.
    InverseLength,          ///< Spring constant inverse proportional to the edge length.
    Harmonic,               ///< Eck's spring constants derived from FEM applied
                            ///< to Laplace equation. The resulting map is the
                            ///< approximate minimizer of the harmonic energy
                            ///< functional. The spring constants can, however,
                            ///< be negative even for convex polygons and thus
                            ///< introduce fold overs in the surface map.
    BarycentricCoordinates, ///< Generalized barycentric coordinates for convex
                            ///< polygons. Requires that the surface is triangulated
                            ///< in such a way that the polygons formed by the
                            ///< vertices adjacent to a given surface node are convex.
                            ///< Coordinates can still be computed if a polygon is
                            ///< concave, but the resulting map may contain fold overs.
    MeanValueCoordinates    ///< Floater's mean value coordinates used to obtain a
                            ///< convex combination map with positive spring constants.
                            ///< These are derived from the Mean Value Theorem of
                            ///< harmonic maps and therefore the resulting map
                            ///< also approximates a harmonic function. However, the
                            ///< approximation is not guaranteed to converge towards
                            ///< the unique continuous solution that minimizes the
                            ///< harmonic energy when the surface is discretized
                            ///< with smaller and smaller triangles. Unlike the
                            ///< map obtained through Eck's or Wachpress' spring
                            ///< constants, the resulting map has no fold overs.
  };

  // ---------------------------------------------------------------------------
  // Attributes

  /// Type of spring constant used to compute surface map
  mirtkPublicAttributeMacro(SpringConstantType, SpringConstant);

  /// Maximum number of iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Tolerance for sparse linear solver
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LinearSpringSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  LinearSpringSurfaceMapper(SpringConstantType = Default);

  /// Copy constructor
  LinearSpringSurfaceMapper(const LinearSpringSurfaceMapper &);

  /// Assignment operator
  LinearSpringSurfaceMapper &operator =(const LinearSpringSurfaceMapper &);

public:

  /// Destructor
  virtual ~LinearSpringSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Remesh surface if necessary
  virtual bool Remesh();

  /// Parameterize surface of input data set
  virtual void Solve();

  /// Whether spring constants are symmetric
  virtual bool Symmetric() const;

  /// Value of spring constant for specified surface edge
  virtual double K(int, int) const;

};


} // namespace mirtk

#endif // MIRTK_LinearSpringSurfaceMapper_H

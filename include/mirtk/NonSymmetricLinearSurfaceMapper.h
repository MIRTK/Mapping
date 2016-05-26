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

#ifndef MIRTK_NonSymmetricLinearSurfaceMapper_H
#define MIRTK_NonSymmetricLinearSurfaceMapper_H

#include "mirtk/LinearSurfaceMapper.h"


namespace mirtk {


/**
 * Obtains surface map as the solution of a non-symmetric system of linear equations
 */
class NonSymmetricLinearSurfaceMapper : public LinearSurfaceMapper
{
  mirtkAbstractMacro(NonSymmetricLinearSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const NonSymmetricLinearSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  NonSymmetricLinearSurfaceMapper();

  /// Copy constructor
  NonSymmetricLinearSurfaceMapper(const NonSymmetricLinearSurfaceMapper &);

  /// Assignment operator
  NonSymmetricLinearSurfaceMapper &operator =(const NonSymmetricLinearSurfaceMapper &);

public:

  /// Destructor
  virtual ~NonSymmetricLinearSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of directed edge (i, j)
  ///
  /// \param[in] i Index of start node.
  /// \param[in] j Index of end node.
  ///
  /// \returns Weight of directed edge (i, j).
  ///
  /// \note This function must be overridden when the base class implementation
  ///       of the Weights function is used. Alternatively, override Weights function.
  virtual double Weight(int i, int j) const;

  /// Weights of edges starting at node i
  ///
  /// \param[in]  i Index of central node.
  /// \param[in]  j Indices of adjacent nodes.
  /// \param[out] w Weights of \p d_i edges (i, j).
  /// \param[in]  d Number of adjacent nodes, i.e., node degree.
  ///
  /// \note The default implementation simply calls the Weight function to compute
  ///       the weight value for each edge independently. When the computation of
  ///       edge weights requires common intermediate computations, it is more
  ///       efficient to compute all edge weights in an overridden subclass
  ///       implementation of this function. The Weight function is then unused.
  virtual void Weights(int i, const int *j, double *w, int d) const;

  /// Construct and solve non-symmetric system of linear equations
  virtual void Solve();

};


} // namespace mirtk

#endif // MIRTK_NonSymmetricLinearSurfaceMapper_H

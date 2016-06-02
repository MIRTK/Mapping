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

#ifndef MIRTK_NearOptimalSurfaceMapper_H
#define MIRTK_NearOptimalSurfaceMapper_H

#include "mirtk/FixedBoundarySurfaceMapper.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Base class of filters that compute a near-optimal intrinsic surface map
 *
 * Surface map filters of this type find the optimal affine combination weights
 * of independently computed conformal and authalic surface maps with given fixed
 * boundary values which minimize a certain non-linear distortion criterion.
 *
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 */
class NearOptimalSurfaceMapper : public FixedBoundarySurfaceMapper
{
  mirtkAbstractMacro(NearOptimalSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of iterations
  ///
  /// When the number of iterations is set to 1 a sparse direct solver is used.
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Tolerance for sparse linear solver
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Computed map values at surface points
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const NearOptimalSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  NearOptimalSurfaceMapper();

  /// Copy constructor
  NearOptimalSurfaceMapper(const NearOptimalSurfaceMapper &);

  /// Assignment operator
  NearOptimalSurfaceMapper &operator =(const NearOptimalSurfaceMapper &);

public:

  /// Destructor
  virtual ~NearOptimalSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Compute map values at interior points
  virtual void ComputeMap();

  /// Assemble output surface map
  virtual void Finalize();

  /// Get component of map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] j Map value component index.
  ///
  /// \return The j-th component of the map value evaluated at the i-th surface point.
  double GetValue(int i, int j = 0) const;

protected:

  /// Compute affine combination weight that minimizes a given distortion measure
  ///
  /// \param[in] u Discrete conformal surface map values.
  /// \param[in] v Discrete authalic  surface map values.
  ///
  /// \returns Affine combination weight \f$\lambda\f$, where final surface map values
  ///          are computed as \f$\lambda u + (1 - \lambda) * v\f$.
  virtual double ComputeLambda(vtkDataArray *u, vtkDataArray *v) const = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double NearOptimalSurfaceMapper::GetValue(int i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}


} // namespace mirtk

#endif // MIRTK_NearOptimalSurfaceMapper_H

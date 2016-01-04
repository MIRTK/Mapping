/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_LinearTetrahedralVolumeParameterizer_H
#define MIRTK_LinearTetrahedralVolumeParameterizer_H

#include <mirtkTetrahedralVolumeParameterizer.h>

#include <vector>


namespace mirtk {


// Forward declarations
class Matrix3x3;


/**
 * Base class of volume re-parameterization filters based on discrete linear operators
 *
 * The reparameterization is the solution of a system of linear equations.
 * The output domain of the volumetric map is currently restricted to 3D.
 * This should not be a practical limitation because these filters are mainly
 * used to re-parameterize a 3D volume, i.e., a map from a 3D domain to another
 * 3D domain.
 */
class LinearTetrahedralVolumeParameterizer : public TetrahedralVolumeParameterizer
{
  mirtkAbstractMacro(LinearTetrahedralVolumeParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of linear solver iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Linear solver tolerance
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Relaxation factor
  mirtkPublicAttributeMacro(double, RelaxationFactor);

  /// Original ID of n-th interior point
  mirtkReadOnlyAttributeMacro(std::vector<int>, InteriorPointId);

  /// Variable/linear equation offset of n-th interior point
  mirtkReadOnlyAttributeMacro(std::vector<int>, InteriorPointPos);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LinearTetrahedralVolumeParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  LinearTetrahedralVolumeParameterizer();

  /// Copy constructor
  LinearTetrahedralVolumeParameterizer(const LinearTetrahedralVolumeParameterizer &);

  /// Assignment operator
  LinearTetrahedralVolumeParameterizer &operator =(const LinearTetrahedralVolumeParameterizer &);

public:

  /// Destructor
  virtual ~LinearTetrahedralVolumeParameterizer();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Parameterize interior points
  virtual void Parameterize();

  /// Solve linear system with operator weights computed using the passed object
  void Solve(const LinearTetrahedralVolumeParameterizer *);

  // ---------------------------------------------------------------------------
  // Auxiliary functions

public:

  /// Calculate operator weight for given tetrahadron
  ///
  /// \param[in] cellId ID of tetrahedron.
  /// \param[in] v0     First  vertex/point of tetrahedron.
  /// \param[in] v1     Second vertex/point of tetrahedron.
  /// \param[in] v2     Third  vertex/point of tetrahedron.
  /// \param[in] v3     Fourth vertex/point of tetrahedron.
  /// \param[in] volume Volume of tetrahedron.
  ///
  /// \return Operator weight contribution of tetrahedron.
  virtual Matrix3x3 GetWeight(vtkIdType cellId,
                              const double v0[3],
                              const double v1[3],
                              const double v2[3],
                              const double v3[3],
                              double       volume) const = 0;

};


} // namespace mirtk

#endif // MIRTK_LinearTetrahedralVolumeParameterizer_H

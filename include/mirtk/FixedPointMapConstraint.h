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

#ifndef MIRTK_FixedPointMapConstraint_H
#define MIRTK_FixedPointMapConstraint_H

#include "mirtk/MapConstraint.h"


namespace mirtk {


/**
 * Assign fixed value to discrete map domain point
 *
 * This map constraint assigns a fixed value to one or more discrete points
 * of the discretized map domain. It is in particular used to assign fixed values
 * to boundary points of a surface mesh when computing a surface map
 * (see BoundaryMapper and FixedBoundarySurfaceMapper).
 */
class NodeValueMapConstraint : public MapConstraint
{
  mirtkObjectMacro(NodeValueMapConstraint);

  // ---------------------------------------------------------------------------
  // Types

  /// Type of surface point ID to boundary point index map
  typedef UnorderedMap<int, int> PointIdToIndexMap;

  // ---------------------------------------------------------------------------
  // Attributes

  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataSet>, Domain);

  /// Number of (boundary) map components
  mirtkPublicAttributeMacro(int, NumberOfComponents);

  /// Fixed map values
  mirtkAttributeMacro(Matrix, FixedValues);
  
  /// Map from fixed point index to surface point ID
  mirtkAttributeMacro(Array<int>, FixedPointIds);

  /// Map from surface point ID to fixed point index
  mirtkAttributeMacro(PointIdToIndexMap, FixedPointIndex);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const FixedBoundarySurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  FixedBoundarySurfaceMapper();

  /// Copy constructor
  FixedBoundarySurfaceMapper(const FixedBoundarySurfaceMapper &);

  /// Assignment operator
  FixedBoundarySurfaceMapper &operator =(const FixedBoundarySurfaceMapper &);

public:

  /// Destructor
  virtual ~FixedBoundarySurfaceMapper();

  // ---------------------------------------------------------------------------
  // Constraints

  /// Reserve memory for fixed constraints
  ///
  /// \param[in] n Number of fixed values.
  void Reserve(int n);

  /// Number of points with a fixed map value
  int NumberOfFixedPoints() const;

  /// Set fixed map value at (boundary) point
  ///
  /// This function adds a Dirichlet surface point constraint.
  ///
  /// \param[in] ptId  Surface point ID.
  /// \param[in] value Fixed map value. NumberOfComponents must be set.
  void SetFixedValue(int ptId, const double *value);

  /// Set fixed map value at (boundary) point
  ///
  /// This function adds a Dirichlet surface point constraint.
  ///
  /// \param[in] ptId  Surface point ID.
  /// \param[in] value Fixed map value.
  void SetFixedValue(int ptId, const Vector &value);

  /// Set fixed map values at (boundary) points
  ///
  /// This function adds multiple Dirichlet surface point constraints at once.
  /// It has the same effect as a call of Reserve followed by n calls of the
  /// FixedValue function to set the constraint values for each fixed point.
  /// The output of a BoundaryMapper can be passed using this function.
  ///
  /// \param[in] ptIds  Surface point IDs.
  /// \param[in] values Fixed values at surface points. Each column of the matrix
  ///                   must contain the map value of the i-th fixed surface point.
  void SetFixedValues(const Array<int> &ptIds, const Matrix &values);

  /// Release unused reserved memory
  virtual void Squeeze();

  /// Clear constraints
  virtual void Clear();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();


};

////////////////////////////////////////////////////////////////////////////////
// Inline
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int FixedBoundarySurfaceMapper::NumberOfFixedPoints() const
{
  return static_cast<int>(_FixedPointIds.size());
}


} // namespace mirtk

#endif // MIRTK_FixedBoundarySurfaceMapper_H

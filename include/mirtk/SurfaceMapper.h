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

#ifndef MIRTK_SurfaceMapper_H
#define MIRTK_SurfaceMapper_H

#include "mirtk/Object.h"

#include "mirtk/Mapping.h"
#include "mirtk/Memory.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"


namespace mirtk {


/**
 * Base class of solvers for the computation of a surface map
 *
 * Solvers of this type compute a map value for each point on the surface
 * of the input point set. The properties of the output map depend on the
 * boundary conditions and the specific solver used to compute the surface map.
 * Examples of surface maps are the assignment of 2D texture coordinates and
 * a bijective mapping of a cortical surface mesh to a disk, square, or sphere,
 * respectively.
 */
class SurfaceMapper : public Object
{
  mirtkAbstractMacro(SurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input surface mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Domain);

  /// Input map values of (fixed) surface points, i.e., boundary conditions
  ///
  /// The entries at interior points of the surface of the input point set are
  /// used to initialize the output map, while entries at the surface boundary
  /// points of the input point set are used as boundary conditions.
  ///
  /// When a \c Mask is given, points with a non-zero mask value are considered
  /// to be either on the boundary or interior points with fixed map values.
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Input);

  /// Input mask of free (mask=0) and fixed (mask!=0) map values
  ///
  /// When no input mask is provided, a mask is created with non-zero value for
  /// all end points of edges with less than two adjacent faces, i.e.,
  /// the boundary edges of the surface mesh.
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Mask);

  /// (Remeshed) surface
  mirtkAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// Mask of free (mask=0) and fixed (mask!=0) map values
  ///
  /// When no input mask is provided, a mask is created with non-zero value for
  /// all end points of edges with less than two adjacent faces, i.e.,
  /// the boundary edges of the surface mesh.
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Fixed);

  /// Map values at surface points
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// IDs of surface points with free map values
  mirtkAttributeMacro(vtkSmartPointer<vtkIdList>, FreePoints);

  /// IDs of surface points with fixed map values
  mirtkAttributeMacro(vtkSmartPointer<vtkIdList>, FixedPoints);

  /// Index of point in set of points with free (index >= 0) or fixed (index < 0) map values
  mirtkAttributeMacro(vtkSmartPointer<vtkIdList>, PointIndex);

  /// Output surface map
  ///
  /// \note The output map is uninitialized! Mapping::Initialize must be
  ///       called before this map can be evaluated at map domain points.
  mirtkReadOnlyAttributeMacro(SharedPtr<Mapping>, Output);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  SurfaceMapper();

  /// Copy constructor
  SurfaceMapper(const SurfaceMapper &);

  /// Assignment operator
  SurfaceMapper &operator =(const SurfaceMapper &);

public:

  /// Destructor
  virtual ~SurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute surface map
  void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Remesh surface if necessary
  virtual bool Remesh();

  /// Compute map values at free surface points
  virtual void Solve() = 0;

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Auxiliaries

protected:

  /// Number of surface points
  int NumberOfPoints() const;

  /// Number of surface points with free map value
  int NumberOfFreePoints() const;

  /// Number of surface points with fixed map value
  int NumberOfFixedPoints() const;

  /// Dimension of codomain of surface map
  int NumberOfComponents() const;

  /// Whether the map value of the specified surface point is fixed
  template <typename IdType>
  bool IsFixedPoint(IdType i) const;

  /// Get index of i-th point with free map value or -1 if map value of point is fixed
  ///
  /// \param[in] i Surface point index.
  ///
  /// \return Index of point in set of points with free map value.
  template <typename IdType>
  int FreePointIndex(IdType i) const;

  /// Get surface point ID of i-th point with free map value
  ///
  /// \param[in] i Index of point in set of points with free map value.
  ///
  /// \return Surface point index.
  int FreePointId(int i) const;

  /// Get index of i-th point with fixed map value or -1 if map value of point is free
  ///
  /// \param[in] i Surface point index.
  ///
  /// \return Index of point in set of points with fixed map value.
  template <typename IdType>
  int FixedPointIndex(IdType i) const;

  /// Get surface point ID of i-th point with fixed map value
  ///
  /// \param[in] i Index of point in set of points with fixed map value.
  ///
  /// \return Surface point index.
  int FixedPointId(int i) const;

  /// Set scalar map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] v Map value.
  template <typename IdType>
  void SetValue(IdType i, double v);

  /// Set component of map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] j Map component index.
  /// \param[in] v Map component value.
  template <typename IdType>
  void SetValue(IdType i, int j, double v);

  /// Get component of map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] j Map value component index.
  ///
  /// \return The j-th component of the map value evaluated at the i-th surface point.
  template <typename IdType>
  double GetValue(IdType i, int j = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfPoints() const
{
  return static_cast<int>(_Surface->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfFixedPoints() const
{
  return static_cast<int>(_FixedPoints->GetNumberOfIds());
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfFreePoints() const
{
  return NumberOfPoints() - NumberOfFixedPoints();
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline bool SurfaceMapper::IsFixedPoint(IdType ptId) const
{
  return _Fixed->GetComponent(static_cast<vtkIdType>(ptId), 0) != .0;
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline int SurfaceMapper::FreePointIndex(IdType ptId) const
{
  const int i = static_cast<int>(_PointIndex->GetId(static_cast<vtkIdType>(ptId)));
  return (i < 0 ? -1 : i);
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::FreePointId(int i) const
{
  return static_cast<int>(_FreePoints->GetId(static_cast<vtkIdType>(i)));
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline int SurfaceMapper::FixedPointIndex(IdType ptId) const
{
  const int i = static_cast<int>(_PointIndex->GetId(static_cast<vtkIdType>(ptId)));
  return (i < 0 ? (-i) - 1 : -1);
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::FixedPointId(int i) const
{
  return static_cast<int>(_FixedPoints->GetId(static_cast<vtkIdType>(i)));
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfComponents() const
{
  return static_cast<int>(_Values->GetNumberOfComponents());
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline void SurfaceMapper::SetValue(IdType i, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), 0, v);
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline void SurfaceMapper::SetValue(IdType i, int j, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), j, v);
}

// -----------------------------------------------------------------------------
template <typename IdType>
inline double SurfaceMapper::GetValue(IdType i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}


} // namespace mirtk

#endif // MIRTK_SurfaceMapper_H

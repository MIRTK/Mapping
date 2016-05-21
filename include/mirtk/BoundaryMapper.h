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

#ifndef MIRTK_BoundaryMapper_H
#define MIRTK_BoundaryMapper_H

#include "mirtk/Object.h"

#include "mirtk/Memory.h"
#include "mirtk/Array.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/Point.h"

#include "mirtk/BoundaryParameterizer.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"


namespace mirtk {


/**
 * Base class of objects which assign map values to a surface boundary
 *
 * Objects of this type take a surface mesh as input and compute a map value
 * for each point on the surface boundary. Such piecewise linear boundary map
 * can be used as boundary condition for a surface map.
 *
 * \sa SurfaceMapper
 */
class BoundaryMapper : public Object
{
  mirtkAbstractMacro(BoundaryMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input surface mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// Edge table of input surface mesh
  mirtkPublicAttributeMacro(SharedPtr<mirtk::EdgeTable>, EdgeTable);

  /// List of user selected points
  ///
  /// How many points need to be selected depends on the specific parameterizer.
  /// In general, no points have to be selected, but selecting boundary points
  /// can be used to orient the curve, i.e., to decide in which direction along
  /// the closed curve the t values are increasing. The curve orientation is
  /// generally defined by the order of the selected points. The first point in
  /// the selection usually defines the boundary point corresponding to the curve
  /// parameter t=0.
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkIdList>, Selection);

  /// Boundary map values
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Closed surface boundary segment parameterizer
  mirtkPublicAttributeMacro(SharedPtr<BoundaryParameterizer>, Parameterizer);

  /// Closed line strips of surface boundary segments
  vtkSmartPointer<vtkCellArray> _BoundarySegments;

  /// Boundary point index of surface point or -1 if point is not on the boundary
  vtkSmartPointer<vtkIdList> _BoundaryPointIndex;

  /// IDs of surface boundary points
  vtkSmartPointer<vtkIdList> _BoundaryPointIds;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundaryMapper &);

  // ---------------------------------------------------------------------------
  // Auxiliary functions -- Initialize function must be called first

public:

  /// Number of map value components
  virtual int NumberOfComponents() const;

  /// Number of surface points
  int NumberOfPoints() const;
 
  /// Number of boundary segments
  int NumberOfBoundarySegments() const;

  /// Number of surface boundary points
  int NumberOfBoundaryPoints() const;
 
  /// Number of points on boundary segment
  int NumberOfBoundaryPoints(int) const;
 
  /// Number of interior surface points
  int NumberOfInteriorPoints() const;

  /// Get indices of surface boundary points making up a boundary segment
  ///
  /// A boundary segment is a closed line strip. The last boundary point
  /// in the returned list of boundary point indices is connected to the
  /// first point. 
  ///
  /// \param[in]  n     Index of boundary segment.
  /// \param[out] ptIds Indices of surface boundary points.
  void BoundaryPointIds(int n, vtkIdList *ptIds) const;

  /// Get indices of surface boundary points making up a boundary segment
  ///
  /// A boundary segment is a closed line strip. The last boundary point
  /// in the returned list of boundary point indices is connected to the
  /// first point. 
  ///
  /// \param[in]  n Index of boundary segment.
  ///
  /// \returns Indices of surface boundary points.
  vtkSmartPointer<vtkIdList> BoundaryPointIds(int n) const;

  /// Get indices of boundary points making up a boundary segment
  ///
  /// A boundary segment is a closed line strip. The last boundary point
  /// in the returned list of boundary point indices is connected to the
  /// first point. 
  ///
  /// \param[in]  n Index of boundary segment.
  /// \param[out] i Indices of boundary points.
  void BoundaryPointIndices(int n, Array<int> &i) const;

  /// Get indices of boundary points making up a boundary segment
  ///
  /// A boundary segment is a closed line strip. The last boundary point
  /// in the returned list of boundary point indices is connected to the
  /// first point. 
  ///
  /// \param[in]  n Index of boundary segment.
  ///
  /// \returns Indices of boundary points.
  Array<int> BoundaryPointIndices(int n) const;

  /// Get index of boundary segment that a surface point belongs to
  ///
  /// \param[in] ptId Index of surface point.
  ///
  /// \returns Index of boundary segment or -1 when point is not a boundary point.
  template <class IdType>
  int BoundarySegmentIndex(IdType ptId) const;

  /// Surface point index of boundary point
  ///
  /// \param[in] i Boundary point index.
  ///
  /// \returns Surface point index of i-th boundary point.
  template <class IndexType>
  int BoundaryPointId(IndexType i) const;
 
  /// Coordinates of boundary point
  ///
  /// \param[in]  i Boundary point index.
  /// \param[out] p Point coordinates.
  template <class IndexType>
  void BoundaryPoint(IndexType i, double p[3]) const;

  /// Coordinates of boundary point
  ///
  /// \param[in]  i Boundary point index.
  ///
  /// \return Point coordinates.
  template <class IndexType>
  Point BoundaryPoint(IndexType i) const;

  /// Index of surface boundary point
  ///
  /// \param[in] ptId Surface point index.
  ///
  /// \returns Index of surface boundary point or -1 if point is not on the boundary.
  template <class IdType>
  int BoundaryPointIndex(IdType ptId) const;
 
  /// Check if surface point is on the boundary
  ///
  /// \param[in] ptId Surface point index.
  ///
  /// \returns Whether a given surface point is on the boundary.
  template <class IdType>
  bool IsBoundaryPoint(IdType ptId) const;

  /// Check if surface point is not on the boundary
  ///
  /// \param[in] ptId Surface point index.
  ///
  /// \returns Whether a given surface point is not on the boundary.
  template <class IdType>
  bool IsInteriorPoint(IdType ptId) const;

  /// Get map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] j Index of map value component.
  ///
  /// \returns Map value component at surface boundary point.
  template <class IndexType>
  double GetBoundaryValue(IndexType i, int j = 0) const;

  /// Set map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] v Value of j-th map component.
  ///
  /// \returns Map value component at surface point.
  template <class IndexType>
  void SetBoundaryValue(IndexType i, double v) const;

  /// Set map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] j Index of map value component.
  /// \param[in] v Value of j-th map component.
  ///
  /// \returns Map value component at surface point.
  template <class IndexType>
  void SetBoundaryValue(IndexType i, int j, double v) const;

  /// Get map value at surface point
  ///
  /// \param[in] ptId Surface point index.
  /// \param[in] j    Index of map value component.
  ///
  /// \returns Map value component at surface point.
  template <class IdType>
  double GetSurfaceValue(IdType ptId, int j = 0) const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  BoundaryMapper();

  /// Copy constructor
  BoundaryMapper(const BoundaryMapper &);

  /// Assignment operator
  BoundaryMapper &operator =(const BoundaryMapper &);

public:

  /// Destructor
  virtual ~BoundaryMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Assign map values to boundary of input surface
  ///
  /// When only a single boundary segment should be mapped or this boundary
  /// mapper processes all boundary segments at once, this function can be
  /// used as a short-hand for calling Initialize, MapBoundary, and Finalize.
  ///
  /// When multiple boundary segments should be processed, call Initialize once
  /// and then repeatedly set the anchor points and change other map attributes
  /// followed by MapBoundary to process the selected boundary segment.
  /// When all boundary segments are processed, generate the surface map by
  /// calling the Finalize function.
  void Run();

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Process boundary segment
  ///
  /// \param[in] n Index of boundary segment.
  virtual void MapBoundary(int n);

  /// Finalize filter execution
  virtual void Finalize();

protected:

  /// Map boundary segment
  ///
  /// \param[in] n         Index of boundary segment.
  /// \param[in] indices   Indices of boundary points forming a closed line strip
  ///                      that discretizes the current surface boundary segment.
  /// \param[in] tvalues   Curve parameter in [0, 1) for each boundary segment point,
  ///                      where the first point has value t=0 and the parameter value
  ///                      for consecutive points is proportional to the distance of
  ///                      the point from the first point along the boundary curve.
  /// \param[in] selection Indices in \p i and \p t arrays corresponding to
  ///                      selected boundary points.
  virtual void MapBoundarySegment(int n, const Array<int>    &indices,
                                         const Array<double> &tvalues,
                                         const Array<int>    &selection) = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int BoundaryMapper::NumberOfPoints() const
{
  return static_cast<int>(_Surface->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int BoundaryMapper::NumberOfBoundaryPoints() const
{
  return static_cast<int>(_BoundaryPointIds->GetNumberOfIds());
}

// -----------------------------------------------------------------------------
inline int BoundaryMapper::NumberOfBoundaryPoints(int n) const
{
  vtkIdType npts, *pts;
  _BoundarySegments->GetCell(n, npts, pts);
  return static_cast<int>(npts);
}

// -----------------------------------------------------------------------------
inline int BoundaryMapper::NumberOfInteriorPoints() const
{
  return NumberOfPoints() - NumberOfBoundaryPoints();
}

// -----------------------------------------------------------------------------
inline void BoundaryMapper::BoundaryPointIds(int n, vtkIdList *ptIds) const
{
  _BoundarySegments->GetCell(n, ptIds);
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkIdList> BoundaryMapper::BoundaryPointIds(int n) const
{
  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  _BoundarySegments->GetCell(n, ptIds);
  return ptIds;
}

// -----------------------------------------------------------------------------
inline Array<int> BoundaryMapper::BoundaryPointIndices(int n) const
{
  Array<int> i;
  BoundaryPointIndices(n, i);
  return i;
}

// -----------------------------------------------------------------------------
template <class IdType>
inline int BoundaryMapper::BoundarySegmentIndex(IdType ptId) const
{
  int i = -1;
  if (_BoundarySegments) {
    vtkSmartPointer<vtkIdList> ptIds;
    for (vtkIdType loc = 0; loc < _BoundarySegments->GetNumberOfCells(); ++loc) {
      _BoundarySegments->GetCell(loc, ptIds);
      if (ptIds->IsId(static_cast<vtkIdType>(ptId)) != -1) {
        i = static_cast<int>(loc);
        break;
      }
    }
  }
  return i;
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline int BoundaryMapper::BoundaryPointId(IndexType i) const
{
  return static_cast<int>(_BoundaryPointIds->GetId(static_cast<vtkIdType>(i)));
}

// -----------------------------------------------------------------------------
template <class IdType>
inline int BoundaryMapper::BoundaryPointIndex(IdType ptId) const
{
  int i = static_cast<int>(_BoundaryPointIndex->GetId(static_cast<vtkIdType>(ptId)));
  return (i < 0 ? -1 : i);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline void BoundaryMapper::BoundaryPoint(IndexType i, double p[3]) const
{
  _Surface->GetPoint(static_cast<vtkIdType>(BoundaryPointId(i)), p);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline Point BoundaryMapper::BoundaryPoint(IndexType i) const
{
  double p[3];
  BoundaryPoint(i, p);
  return Point(p);
}

// -----------------------------------------------------------------------------
template <class IdType>
inline bool BoundaryMapper::IsBoundaryPoint(IdType ptId) const
{
  return BoundaryPointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
template <class IdType>
inline bool BoundaryMapper::IsInteriorPoint(IdType ptId) const
{
  return !IsBoundaryPoint(ptId);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline double BoundaryMapper::GetBoundaryValue(IndexType i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline void BoundaryMapper::SetBoundaryValue(IndexType i, double v) const
{
  return _Values->SetComponent(static_cast<vtkIdType>(i), 0, v);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline void BoundaryMapper::SetBoundaryValue(IndexType i, int j, double v) const
{
  return _Values->SetComponent(static_cast<vtkIdType>(i), j, v);
}

// -----------------------------------------------------------------------------
template <class IdType>
inline double BoundaryMapper::GetSurfaceValue(IdType ptId, int j) const
{
  const int i = BoundaryPointIndex(ptId);
  if (i < 0) {
    if (_Values) {
      return _Values->GetComponent(static_cast<vtkIdType>(ptId), j);
    } else {
      return .0;
    }
  }
  return GetBoundaryValue(i, j);
}


} // namespace mirtk

#endif // MIRTK_BoundaryMapper_H

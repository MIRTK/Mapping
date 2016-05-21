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

#ifndef MIRTK_BoundaryParameterizer_H
#define MIRTK_BoundaryParameterizer_H

#include "mirtk/Object.h"
#include "mirtk/Array.h"
#include "mirtk/Math.h"
#include "mirtk/Vector.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkIdList.h"


namespace mirtk {


/**
 * Base class of closed boundary curve parameterizers
 *
 * A boundary parameterizer is used to assign each point of a closed surface
 * boundary segment a parameter value t in [0, 1), with strictly increasing
 * (or decreasing) values along the boundary curve. This parameterization is
 * used by a boundary mapper to assign map values to the boundary points.
 * For example, the BoundaryToSquareMapper maps points with a t value in
 * [0, .25) to the first side of the square, points with a t value in [.25, .5)
 * to the second side, and so on. Note that the t value is proportional to the
 * distance of the point from the curve point at t=0 along the mapped curve.
 * Note futher that the boundary point with index 0 need not be the boundary
 * point with parameter value t=0. It is also not required that there is a
 * discrete curve point with parameter t=0.
 *
 * Specialized boundary parameterizers may assign parameter values to the
 * boundary segment of a given surface mesh not only based on this surface
 * and possibly related imaging data, but also other surface meshes with
 * known or to be determined correspondences between surface boundaries.
 * Such parameterizer will assign equal t values to corresponding boundary
 * points such that these points are mapped to the same value by the boundary
 * mapper. Other parameterizers may take an ordered/labeled list of manually
 * selected boundary points as input and assign parameter values based on
 * this selection.
 *
 * In a nutshell, a boundary parameterizer indirectly establishes
 * correspondences between boundary points of 2D manifolds, while the boundary
 * mapper itself defines the shape of the map codomain.
 *
 * \sa BoundaryMapper
 */
class BoundaryParameterizer : public Object
{
  mirtkAbstractMacro(BoundaryParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Ordered list of indices of points making up the surface boundary segment
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkIdList>, Boundary);

  /// Coordinates of surface points
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPoints>, Points);

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

  /// Selected points of surface boundary segment
  mirtkAttributeMacro(vtkSmartPointer<vtkIdList>, BoundarySelection);

  /// Parameter value at the boundary points
  mirtkReadOnlyAttributeMacro(Array<double>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundaryParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  BoundaryParameterizer();

  /// Copy constructor
  BoundaryParameterizer(const BoundaryParameterizer &);

  /// Assignment operator
  BoundaryParameterizer &operator =(const BoundaryParameterizer &);

public:

  /// Destructor
  virtual ~BoundaryParameterizer();

  /// New copy of this parameterizer
  virtual BoundaryParameterizer *NewCopy() const = 0;

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Number of boundary curve points
  int NumberOfBoundaryPoints() const;

  /// Number of selected curve points
  int NumberOfSelectedPoints() const;

  /// Convert any index to a boundary point index in [0, n)
  template <class IndexType>
  IndexType ToBoundaryPointIndexRange(IndexType i) const;

  /// Get surface point coordinates of i-th boundary segment point
  ///
  /// \param[in]  i Index of boundary segment point.
  ///               The index is taken modulo the number of boundary
  ///               points and negative values count from the end of
  ///               the list of boundary points.
  /// \param[out] p Boundary point coordinates.
  template <class IndexType>
  void BoundaryPoint(IndexType i, double p[3]) const;
 
  /// Get surface point ID of i-th boundary segment point
  ///
  /// \param[in] i Index of boundary segment point.
  ///              The index is taken modulo the number of boundary
  ///              points and negative values count from the end of
  ///              the list of boundary points.
  ///
  /// \return Index of corresponding surface point.
  template <class IndexType>
  int BoundaryPointId(IndexType i) const;

  /// Get boundary segment point index
  ///
  /// \param[in] ptId Index of surface point.
  ///
  /// \return Index of corresponding boundary segment point.
  /// \retval -1 If point is not on this boundary segment.
  template <class IdType>
  int BoundaryPointIndex(IdType ptId) const;

  /// Check if surface point is on the boundary segment
  ///
  /// \param[in] ptId Index of surface point.
  ///
  /// \return Whether surface point is on this boundary segment.
  template <class IdType>
  bool IsBoundaryPoint(IdType ptId) const;

  /// Get surface point ID of i-th selected boundary point
  ///
  /// \param[in] i Index of selected point.
  ///
  /// \return Index of corresponding surface point.
  template <class IndexType>
  int SelectedPointId(IndexType i) const;

  /// Get boundary segment point index of i-th selected boundary point
  ///
  /// \param[in] ptId Index of selected point.
  ///
  /// \return Index of corresponding boundary segment point.
  template <class IndexType>
  int SelectedPointIndex(IndexType ptId) const;

  /// Check if boundary segment point is selected
  ///
  /// \param[in] i Index of boundary segment point.
  ///              The index is taken modulo the number of boundary
  ///              points and negative values count from the end of
  ///              the list of boundary points.
  ///
  /// \return Whether boundary segment point is selected.
  template <class IndexType>
  bool IsSelected(IndexType i) const;

protected:

  /// Distance between points
  static double Distance(double p1[3], double p2[3]);

  /// Lengths of boundary segment edges
  Vector BoundaryEdgeLengths() const;

  /// Length of next boundary edge
  ///
  /// \param[in] i Index of boundary segment point.
  /// \param[in] l Precomputed edge lenghts.
  ///
  template <class IndexType>
  double BoundaryEdgeLength(IndexType i, const Vector *l = nullptr) const;

  /// Length of next boundary edge
  ///
  /// \param[in] i  Index of boundary segment point.
  /// \param[in] di Index increment +1 or -1, i.e., traversal direction.
  /// \param[in] l  Precomputed edge lenghts.
  ///
  template <class IndexType>
  double BoundaryEdgeLength(IndexType i, IndexType di, const Vector *l = nullptr) const;

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Parameterize boundary segment
  virtual void Run();

protected:

  /// Initialize parameterizer after input and parameters are set
  virtual void Initialize();

  /// Parameterize boundary curve
  virtual void Parameterize() = 0;

  /// Finalize parameterization
  virtual void Finalize();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int BoundaryParameterizer::NumberOfBoundaryPoints() const
{
  return static_cast<int>(_Boundary->GetNumberOfIds());
}

// -----------------------------------------------------------------------------
inline int BoundaryParameterizer::NumberOfSelectedPoints() const
{
  return static_cast<int>(_BoundarySelection->GetNumberOfIds());
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline IndexType BoundaryParameterizer::ToBoundaryPointIndexRange(IndexType i) const
{
  const IndexType max = static_cast<IndexType>(NumberOfBoundaryPoints() - 1);
  if (i < 0) {
    i = -i;
    int n = i / max;
    int m = i - n * max;
    if (n & 1) i = max - m;
    else       i = m;
  } else if (i > max) {
    i -= max;
    int n = i / max;
    int m = i - n * max;
    if (n & 1) i = m;
    else       i = max - m;
  }
  return i;
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline void BoundaryParameterizer::BoundaryPoint(IndexType i, double p[3]) const
{
  vtkIdType pos = static_cast<vtkIdType>(ToBoundaryPointIndexRange(i));
  _Points->GetPoint(_Boundary->GetId(pos), p);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline int BoundaryParameterizer::BoundaryPointId(IndexType i) const
{
  vtkIdType pos = static_cast<vtkIdType>(ToBoundaryPointIndexRange(i));
  return static_cast<int>(_Boundary->GetId(pos));
}

// -----------------------------------------------------------------------------
template <class IdType>
inline int BoundaryParameterizer::BoundaryPointIndex(IdType ptId) const
{
  int i = static_cast<int>(_Boundary->IsId(static_cast<vtkIdType>(ptId)));
  return i;
}

// -----------------------------------------------------------------------------
template <class IdType>
inline bool BoundaryParameterizer::IsBoundaryPoint(IdType ptId) const
{
  return BoundaryPointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline int BoundaryParameterizer::SelectedPointId(IndexType i) const
{
  return static_cast<int>(_BoundarySelection->GetId(static_cast<vtkIdType>(i)));
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline int BoundaryParameterizer::SelectedPointIndex(IndexType i) const
{
  return BoundaryPointIndex(SelectedPointId(i));
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline bool BoundaryParameterizer::IsSelected(IndexType i) const
{
  vtkIdType ptId = static_cast<vtkIdType>(BoundaryPointId(i));
  return _BoundarySelection->IsId(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline double BoundaryParameterizer::Distance(double p1[3], double p2[3])
{
  double dx = p2[0] - p1[0];
  double dy = p2[1] - p1[1];
  double dz = p2[2] - p1[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline double BoundaryParameterizer::BoundaryEdgeLength(IndexType i, IndexType di, const Vector *l) const
{
  if (l == nullptr) {
    double p1[3], p2[3];
    BoundaryPoint(i,      p1);
    BoundaryPoint(i + di, p2);
    return Distance(p1, p2);
  } else {
    return l(ToBoundaryPointIndexRange(di < 0 ? i-1 : i));
  }
}

// -----------------------------------------------------------------------------
template <class IndexType>
inline double BoundaryParameterizer::BoundaryEdgeLength(IndexType i, const Vector *l) const
{
  return BoundaryEdgeLength(i, +1, l);
}


} // namespace mirtk

#endif // MIRTK_BoundaryParameterizer_H

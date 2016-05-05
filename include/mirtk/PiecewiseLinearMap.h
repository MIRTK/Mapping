/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_PiecewiseLinearMap_H
#define MIRTK_PiecewiseLinearMap_H

#include "mirtk/Mapping.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkDataArray.h"
#include "vtkAbstractCellLocator.h"


namespace mirtk {


/**
 * Piecewise linear map defined at mesh nodes
 *
 * This map is defined by its values at discrete mesh nodes.
 * Intermediate values are interpolated using the weights of the respective
 * mesh cell which this point belongs to. If the point is not contained in
 * any mesh cell, the point is mapped to a constant outside value.
 *
 * A surface map is commonly represented by the texture coordinates of a
 * triangular surface mesh, while a tetrahedral mesh is commonly used to
 * parameterize a volumetric map. These maps are usually computed using a
 * finite element method (FEM).
 */
class PiecewiseLinearMap : public Mapping
{
  mirtkObjectMacro(PiecewiseLinearMap);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Mesh which discretizes the domain of this piecewise linear map
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataSet>, Domain);

  /// Point data array with map values for each mesh point
  ///
  /// If not set, the active TCOORDS, VECTORS, or SCALARS array is used.
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Locates cell within which a given point lies
  mirtkAttributeMacro(vtkSmartPointer<vtkAbstractCellLocator>, Locator);

  /// Maximum number cell points
  mirtkAttributeMacro(int, MaxCellSize);

  /// Squared distance tolerance used to locate cells
  /// \note Unused argument of vtkCellLocator::FindCell (as of VTK <= 7.0).
  static const double _Tolerance2;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const PiecewiseLinearMap &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  PiecewiseLinearMap();

  /// Copy constructor
  PiecewiseLinearMap(const PiecewiseLinearMap &);

  /// Assignment operator
  PiecewiseLinearMap &operator =(const PiecewiseLinearMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual Mapping *NewCopy() const;

  /// Destructor
  virtual ~PiecewiseLinearMap();

  // ---------------------------------------------------------------------------
  // Map domain

  // Import other overloads
  using Mapping::BoundingBox;

  /// Get minimum axes-aligned bounding box of map domain
  ///
  /// \param[out] x1 Lower bound of map domain along x axis.
  /// \param[out] y1 Lower bound of map domain along y axis.
  /// \param[out] z1 Lower bound of map domain along z axis.
  /// \param[out] x2 Upper bound of map domain along x axis.
  /// \param[out] y2 Upper bound of map domain along y axis.
  /// \param[out] z2 Upper bound of map domain along z axis.
  virtual void BoundingBox(double &x1, double &y1, double &z1,
                           double &x2, double &y2, double &z2) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using Mapping::Evaluate;

  /// Dimension of codomain, i.e., number of output values
  virtual int NumberOfComponents() const;

  /// Evaluate map at a given mesh point
  ///
  /// \param[out] v Map value.
  /// \param[in]  i Point index.
  virtual void Evaluate(double *v, int i) const;

  /// Evaluate map at a given mesh point
  ///
  /// \param[in] i Point index.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value at the given mesh point.
  virtual double Evaluate(int i, int l = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  x Coordinate of point along x axis at which to evaluate map.
  /// \param[in]  y Coordinate of point along y axis at which to evaluate map.
  /// \param[in]  z Coordinate of point along z axis at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] x Coordinate of point along x axis at which to evaluate map.
  /// \param[in] y Coordinate of point along y axis at which to evaluate map.
  /// \param[in] z Coordinate of point along z axis at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluate at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Read map from file
  virtual bool Read(const char *);

  /// Write map to file
  virtual bool Write(const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline void PiecewiseLinearMap::Evaluate(double *v, int i) const
{
  _Values->GetTuple(static_cast<vtkIdType>(i), v);
}

// -----------------------------------------------------------------------------
inline double PiecewiseLinearMap::Evaluate(int i, int l) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), l);
}


} // namespace mirtk

#endif // MIRTK_PiecewiseLinearMap_H

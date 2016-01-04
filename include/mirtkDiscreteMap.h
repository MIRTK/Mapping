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

#ifndef MIRTK_DiscreteMap_H
#define MIRTK_DiscreteMap_H

#include <mirtkVolumetricMap.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkDataArray.h>
#include <vtkAbstractCellLocator.h>


namespace mirtk {


/**
 * Volumetric map defined at mesh nodes
 *
 * This volumetric map is defined by its values at discrete mesh nodes.
 * Values between nodes are interpolated using the weights of the respective
 * mesh cell which this point belongs to. If the point is not contained in
 * any mesh cell, the point is mapped to a constant outside value.
 *
 * A tetrahedral mesh is commonly used to parameterize a discrete volumetric map.
 */
class DiscreteMap : public VolumetricMap
{
  mirtkObjectMacro(DiscreteMap);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input mesh which defines this volumetric map
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, Input);

  /// Point data array with map values for each input mesh point
  ///
  /// If NULL, the active SCALARS or VECTORS point data array is used.
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Locates cell which is closest to a given point
  mirtkAttributeMacro(vtkSmartPointer<vtkAbstractCellLocator>, Locator);

  /// Squared distance tolerance used to located cells
  /// \note Currently (VTK <= 6.3) unused argument of vtkCellLocator::FindCell.
  static const double _Tolerance2;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const DiscreteMap &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  DiscreteMap();

  /// Copy constructor
  DiscreteMap(const DiscreteMap &);

  /// Assignment operator
  DiscreteMap &operator =(const DiscreteMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual VolumetricMap *NewCopy() const;

  /// Destructor
  virtual ~DiscreteMap();

  // ---------------------------------------------------------------------------
  // Input domain

  // Import other overloads
  using VolumetricMap::BoundingBox;

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] x1 Lower bound of input domain along x axis.
  /// \param[out] y1 Lower bound of input domain along y axis.
  /// \param[out] z1 Lower bound of input domain along z axis.
  /// \param[out] x2 Upper bound of input domain along x axis.
  /// \param[out] y2 Upper bound of input domain along y axis.
  /// \param[out] z2 Upper bound of input domain along z axis.
  virtual void BoundingBox(double &x1, double &y1, double &z1,
                           double &x2, double &y2, double &z2) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using VolumetricMap::Evaluate;

  /// Dimension of output domain, i.e., number of output values
  virtual int NumberOfComponents() const;

  /// Evaluate vector-valued volumetric map at the given point
  ///
  /// \param[out] v Value of volumetric map.
  /// \param[in]  x Coordinate of point at which to evaluate volumetric map.
  /// \param[in]  y Coordinate of point at which to evaluate volumetric map.
  /// \param[in]  z Coordinate of point at which to evaluate volumetric map.
  ///
  /// \returns Whether input point is within input domain.
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const;

  /// Evaluate scalar volumetric map at the given point
  ///
  /// \param[in] x Coordinate of point at which to evaluate volumetric map.
  /// \param[in] y Coordinate of point at which to evaluate volumetric map.
  /// \param[in] z Coordinate of point at which to evaluate volumetric map.
  /// \param[in] l Index of scalar volumetric map component.
  ///
  /// \returns Scalar value of volumetric map.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Read volumetric map from file
  virtual bool Read(const char *);

  /// Write volumetric map to file
  virtual bool Write(const char *) const;

};


} // namespace mirtk

#endif // MIRTK_DiscreteMap_H

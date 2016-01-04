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

#ifndef MIRTK_VolumetricMap_H
#define MIRTK_VolumetricMap_H

#include <mirtkObject.h>

#include <mirtkPoint.h>
#include <mirtkCfstream.h>
#include <mirtkImageAttributes.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>


namespace mirtk {


// Forward declaration of possible output map type
template <class TVoxel> class GenericImage;


/**
 * Volumetric map
 *
 * A volumetric map assigns each point in a 3D volume a n-D target value.
 */
class VolumetricMap : public Object
{
  mirtkAbstractMacro(VolumetricMap);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Value assigned to points outside the input domain
  mirtkPublicAttributeMacro(double, OutsideValue);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const VolumetricMap &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  VolumetricMap();

  /// Copy constructor
  VolumetricMap(const VolumetricMap &);

  /// Assignment operator
  VolumetricMap &operator =(const VolumetricMap &);

public:

  /// Read volumetric map from file
  static VolumetricMap *New(const char *);

  /// Destructor
  virtual ~VolumetricMap();

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual VolumetricMap *NewCopy() const = 0;

  // ---------------------------------------------------------------------------
  // Input domain

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] x1 Lower bound of input domain along x axis.
  /// \param[out] y1 Lower bound of input domain along y axis.
  /// \param[out] z1 Lower bound of input domain along z axis.
  /// \param[out] x2 Upper bound of input domain along x axis.
  /// \param[out] y2 Upper bound of input domain along y axis.
  /// \param[out] z2 Upper bound of input domain along z axis.
  virtual void BoundingBox(double &x1, double &y1, double &z1,
                           double &x2, double &y2, double &z2) const = 0;

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] bounds Bounds of input domain in VTK order, i.e.,
  ///                    [x1, x2, y1, y2, z1, z2].
  void BoundingBox(double bounds[6]) const;

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] p1 Lower-left-front corner of input domain bounding box.
  /// \param[out] p2 Upper-right-back corner of input domain bounding box.
  void BoundingBox(Point &p1, Point &p2) const;

  /// Get volumetric grid attributes of input domain
  ///
  /// \param[in] nx Number of grid points in x direction.
  /// \param[in] ny Number of grid points in y direction.
  /// \param[in] nz Number of grid points in z direction.
  ImageAttributes Attributes(int nx, int ny = 0, int nz = 0) const;

  /// Get volumetric grid attributes of input domain
  ///
  /// If no grid spacing is specified, the length of the bounding box diagonal
  /// is dividied by 256.
  ///
  /// \param[in] dx Grid spacing along x axis.
  /// \param[in] dy Grid spacing along y axis.
  /// \param[in] dz Grid spacing along z axis.
  ImageAttributes Attributes(double dx = .0, double dy = .0, double dz = .0) const;

  // ---------------------------------------------------------------------------
  // Evaluation

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
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const = 0;

  /// Evaluate vector-valued volumetric map at the given point
  ///
  /// \param[out] v Value of volumetric map.
  /// \param[in]  p Point at which to evaluate volumetric map.
  ///
  /// \returns Whether input point is within input domain.
  bool Evaluate(double *v, const double p[3]) const;

  /// Evaluate vector-valued volumetric map at the given point
  ///
  /// \param[out] v Value of volumetric map.
  /// \param[in]  p Point at which to evaluate volumetric map.
  ///
  /// \returns Whether input point is within input domain.
  bool Evaluate(double *v, const Point &p) const;

  /// Evaluate scalar volumetric map at the given point
  ///
  /// \param[in] x Coordinate of point at which to evaluate volumetric map.
  /// \param[in] y Coordinate of point at which to evaluate volumetric map.
  /// \param[in] z Coordinate of point at which to evaluate volumetric map.
  /// \param[in] l Index of scalar volumetric map component.
  ///
  /// \returns Scalar value of volumetric map.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

  /// Evaluate scalar volumetric map at the given point
  ///
  /// \param[in] p Point at which to evaluate volumetric map.
  /// \param[in] l Index of scalar volumetric map component.
  ///
  /// \returns Scalar value of volumetric map.
  double Evaluate(const double p[3], int l = 0) const;

  /// Evaluate scalar volumetric map at the given point
  ///
  /// \param[in] p Point at which to evaluate volumetric map.
  /// \param[in] l Index of scalar volumetric map component.
  ///
  /// \returns Scalar value of volumetric map.
  double Evaluate(const Point &, int l = 0) const;

  /// Evaluate volumetric map at each point of a regular lattice
  virtual void Evaluate(GenericImage<float> &f, int l = 0, vtkSmartPointer<vtkPointSet> = NULL) const;

  /// Evaluate volumetric map at each point of a regular lattice
  virtual void Evaluate(GenericImage<double> &f, int l = 0, vtkSmartPointer<vtkPointSet> = NULL) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Read volumetric map from file
  virtual bool Read(const char *);

  /// Write volumetric map to file
  virtual bool Write(const char *) const;

protected:

  /// Read attributes of volumetric map from file stream
  virtual void ReadMap(Cifstream &);

  /// Write attributes of volumetric map to file stream
  virtual void WriteMap(Cofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Input domain
// =============================================================================

// -----------------------------------------------------------------------------
inline void VolumetricMap::BoundingBox(double bounds[6]) const
{
  this->BoundingBox(bounds[0], bounds[2], bounds[4],
                    bounds[1], bounds[3], bounds[5]);
}

// -----------------------------------------------------------------------------
inline void VolumetricMap::BoundingBox(Point &p1, Point &p2) const
{
  this->BoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline int VolumetricMap::NumberOfComponents() const
{
  return 1;
}

// -----------------------------------------------------------------------------
inline bool VolumetricMap::Evaluate(double *v, const double p[3]) const
{
  return this->Evaluate(v, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool VolumetricMap::Evaluate(double *v, const Point &p) const
{
  return this->Evaluate(v, p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline double VolumetricMap::Evaluate(double x, double y, double z, int l) const
{
  double * const v = new double[this->NumberOfComponents()];
  this->Evaluate(v, x, y, z);
  const double s = v[l];
  delete[] v;
  return s;
}

// -----------------------------------------------------------------------------
inline double VolumetricMap::Evaluate(const double p[3], int l) const
{
  return this->Evaluate(p[0], p[1], p[2], l);
}

// -----------------------------------------------------------------------------
inline double VolumetricMap::Evaluate(const Point &p, int l) const
{
  return this->Evaluate(p._x, p._y, p._z, l);
}


} // namespace mirtk

#endif // MIRTK_VolumetricMap_H

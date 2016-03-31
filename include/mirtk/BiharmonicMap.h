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

#ifndef MIRTK_BiharmonicMap_H
#define MIRTK_BiharmonicMap_H

#include "mirtk/HarmonicMap.h"


namespace mirtk {


/**
 * Biharmonic volumetric map
 */
class BiharmonicMap : public HarmonicMap
{
  mirtkObjectMacro(BiharmonicMap);

public:

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Biharmonic kernel function
  ///
  /// \param[in] d Distance of point to source point.
  ///
  /// \returns Biarmonic kernel function value.
  static double B(double d);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  BiharmonicMap();

  /// Copy constructor
  BiharmonicMap(const BiharmonicMap &);

  /// Assignment operator
  BiharmonicMap &operator =(const BiharmonicMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual VolumetricMap *NewCopy() const;

  /// Destructor
  virtual ~BiharmonicMap();

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using HarmonicMap::Evaluate;

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

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double BiharmonicMap::B(double d)
{
  return d / (8.0 * Pi());
}


} // namespace mirtk

#endif // MIRTK_BiharmonicMap_H

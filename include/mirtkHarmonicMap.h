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

#ifndef MIRTK_HarmonicMap_H
#define MIRTK_HarmonicMap_H

#include <mirtkFundamentalMap.h>

#include <mirtkPointSet.h>
#include <mirtkVector.h>


namespace mirtk {


/**
 * Harmonic volumetric map
 */
class HarmonicMap : public FundamentalMap
{
  mirtkObjectMacro(HarmonicMap);

public:

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Harmonic kernel function
  ///
  /// \param[in] d Distance of point to source point.
  ///
  /// \returns Harmonic kernel function value.
  static double H(double d);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  HarmonicMap();

  /// Copy constructor
  HarmonicMap(const HarmonicMap &);

  /// Assignment operator
  HarmonicMap &operator =(const HarmonicMap &);

  /// Make deep copy of this volumetric map
  virtual VolumetricMap *NewCopy() const;

  /// Destructor
  virtual ~HarmonicMap();

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using FundamentalMap::Evaluate;

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
inline double HarmonicMap::H(double d)
{
  return .25 / (d * Pi());
}


} // namespace mirtk

#endif // MIRTK_HarmonicMap_H

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

#ifndef MIRTK_HarmonicTetrahedralVolumeParameterizer_H
#define MIRTK_HarmonicTetrahedralVolumeParameterizer_H

#include <mirtkLinearTetrahedralVolumeParameterizer.h>


namespace mirtk {


/**
 * Discrete volumetric harmonic map
 *
 * This implementation is based on
 *
 *   Paillé & Poulin (2012), As-conformal-as-possible discrete volumetric mapping,
 *   Computers and Graphics (Pergamon), 36(5), 427–433. doi:10.1016/j.cag.2012.03.014
 *
 * The discrete volumetric harmonic map was first presented in
 *
 *   Wang et al. (2004), Volumetric harmonic map,
 *   Communications in Information and Systems, 3(3), 191–202.
 */
class HarmonicTetrahedralVolumeParameterizer
: public LinearTetrahedralVolumeParameterizer
{
  mirtkObjectMacro(HarmonicTetrahedralVolumeParameterizer);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  HarmonicTetrahedralVolumeParameterizer();

  /// Copy constructor
  HarmonicTetrahedralVolumeParameterizer(const HarmonicTetrahedralVolumeParameterizer &);

  /// Assignment operator
  HarmonicTetrahedralVolumeParameterizer &operator =(const HarmonicTetrahedralVolumeParameterizer &);

  /// Destructor
  virtual ~HarmonicTetrahedralVolumeParameterizer();

  // ---------------------------------------------------------------------------
  // Auxiliary functions

protected:

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
                              double       volume) const;

};


} // namespace mirtk

#endif // MIRTK_HarmonicTetrahedralVolumeParameterizer_H 

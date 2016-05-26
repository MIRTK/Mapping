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

#ifndef MIRTK_UniformBoundaryParameterizer_H
#define MIRTK_UniformBoundaryParameterizer_H

#include "mirtk/BoundaryParameterizer.h"


namespace mirtk {


/**
 * Boundary curve parameterization with uniform distance of curve points
 *
 * This boundary curve parameterization is referred to as uniform
 * parameterization in Floater (1997).
 *
 * - Floater (1997). Parametrization and smooth approximation of surface triangulations.
 *   Computer Aided Geometric Design, 14(3), 231–250.
 */
class UniformBoundaryParameterizer : public BoundaryParameterizer
{
  mirtkObjectMacro(UniformBoundaryParameterizer);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  UniformBoundaryParameterizer();

  /// Copy constructor
  UniformBoundaryParameterizer(const UniformBoundaryParameterizer &);

  /// Assignment operator
  UniformBoundaryParameterizer &operator =(const UniformBoundaryParameterizer &);

  /// Destructor
  virtual ~UniformBoundaryParameterizer();

  /// New copy of this parameterizer
  virtual BoundaryParameterizer *NewCopy() const;

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Parameterize boundary curve
  virtual void Parameterize();

};


} // namespace mirtk

#endif // MIRTK_UniformBoundaryParameterizer_H

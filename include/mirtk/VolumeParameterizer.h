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

#ifndef MIRTK_VolumeParameterizer_H
#define MIRTK_VolumeParameterizer_H

#include "mirtk/Object.h"

#include "mirtk/Mapping.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Base class of filters which parameterize the interior of a
 * piecewise linear complex (PLC)
 */
class VolumeParameterizer : public Object
{
  mirtkAbstractMacro(VolumeParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input point set (e.g., surface mesh or tetrahedral mesh)
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, InputSet);

  /// Input boundary map
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, InputMap);

  /// Boundary surface of input point set
  mirtkAttributeMacro(vtkSmartPointer<vtkPolyData>, Boundary);

  /// Boundary surface map
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, BoundaryMap);

  /// Volumetric map
  ///
  /// \note The output map is uninitialized! Mapping::Initialize must be
  ///       called before this map can be evaluated at map domain points.
  mirtkReadOnlyComponentMacro(Mapping, OutputMap);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const VolumeParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  VolumeParameterizer();

  /// Copy constructor
  VolumeParameterizer(const VolumeParameterizer &);

  /// Assignment operator
  VolumeParameterizer &operator =(const VolumeParameterizer &);

public:

  /// Destructor
  virtual ~VolumeParameterizer();

  /// Dimension of codomain of volumetric map
  int OutputDimension() const;

  /// Get volumetric map and set _OutputMap to nullptr
  ///
  /// \note The returned object has to be deleted by the caller.
  ///       Use OutputMap() instead to keep ownership with this class.
  Mapping *GetOutputMap();

  // ---------------------------------------------------------------------------
  // Execution

  /// Parameterize interior of input data set
  void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Initialize boundary surface with corresponding boundary map as point data
  virtual void InitializeBoundary(vtkPointSet *, vtkDataArray *);

  /// Parameterize interior of input data set
  virtual void Parameterize() = 0;

  /// Finalize filter execution
  virtual void Finalize();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int VolumeParameterizer::OutputDimension() const
{
  return _InputMap ? static_cast<int>(_InputMap->GetNumberOfComponents()) : 0;
}

// -----------------------------------------------------------------------------
inline Mapping *VolumeParameterizer::GetOutputMap()
{
  Mapping *map = _OutputMap;
  _OutputMap = nullptr;
  return map;
}


} // namespace mirtk

#endif // MIRTK_VolumeParameterizer_H

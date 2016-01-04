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

#ifndef MIRTK_VolumeParameterizer_H
#define MIRTK_VolumeParameterizer_H

#include <mirtkObject.h>

#include <mirtkVolumetricMap.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>


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

  /// Volumetric output map
  /// \note The output map is uninitialized!
  mirtkReadOnlyComponentMacro(VolumetricMap, OutputMap);

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

  /// Dimension of output domain of volumetric map
  int OutputDimension() const;

  /// Get volumetric output map and set _OutputMap to NULL
  ///
  /// \note The returned object has to be deleted by the caller.
  ///       Use OutputMap() instead to keep ownership with this class.
  VolumetricMap *GetOutputMap();

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
inline VolumetricMap *VolumeParameterizer::GetOutputMap()
{
  VolumetricMap *map = _OutputMap;
  _OutputMap = NULL;
  return map;
}


} // namespace mirtk

#endif // MIRTK_VolumeParameterizer_H

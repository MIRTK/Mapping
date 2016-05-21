# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK Mapping library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/${PROJECT_NAME}Export.h
  # Map types
  Mapping.h
  MeshlessMap.h
  MeshlessHarmonicMap.h
  MeshlessBiharmonicMap.h
  PiecewiseLinearMap.h
  # Surface boundary curve parameterization
  BoundaryParameterizer.h
  LinearBoundaryParameterizer.h
  # Surface boundary mapping
  BoundaryMapper.h
  BoundaryToDiskMapper.h
  BoundaryToSquareMapper.h
  BoundaryToPolygonMapper.h
  # Surface mapping
  SurfaceMapper.h
  LinearSpringSurfaceMapper.h
  # Volumetric mapping
  VolumeMapper.h
  TetrahedralMeshMapper.h
  LinearTetrahedralMeshMapper.h
  HarmonicTetrahedralMeshMapper.h
  AsConformalAsPossibleMapper.h
  MeshlessVolumeMapper.h
  MeshlessHarmonicVolumeMapper.h
)

set(SOURCES
  # Map types
  Mapping.cc
  MeshlessMap.cc
  MeshlessHarmonicMap.cc
  MeshlessBiharmonicMap.cc
  PiecewiseLinearMap.cc
  # Surface boundary curve parameterization
  BoundaryParameterizer.cc
  LinearBoundaryParameterizer.cc
  # Surface boundary mapping
  BoundaryMapper.cc
  BoundaryToDiskMapper.cc
  BoundaryToSquareMapper.cc
  BoundaryToPolygonMapper.cc
  # Surface mapping
  SurfaceMapper.cc
  LinearSpringSurfaceMapper.cc
  # Volumetric mapping
  VolumeMapper.cc
  TetrahedralMeshMapper.cc
  LinearTetrahedralMeshMapper.cc
  HarmonicTetrahedralMeshMapper.cc
  AsConformalAsPossibleMapper.cc
  MeshlessVolumeMapper.cc
  MeshlessHarmonicVolumeMapper.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibPointSet
  LibIO
  ${VTK_LIBRARIES}
)

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

mirtk_add_library()

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
  # Mappings
  Mapping.h
  PiecewiseLinearMap.h
  FundamentalMap.h
  HarmonicFundamentalMap.h
  BiharmonicFundamentalMap.h
  # Surface mapping
  # ...
  # Volumetric mapping
  VolumeParameterizer.h
  TetrahedralVolumeParameterizer.h
  LinearTetrahedralVolumeParameterizer.h
  FundamentalVolumeParameterizer.h
  HarmonicTetrahedralVolumeParameterizer.h
  HarmonicFundamentalVolumeParameterizer.h
  AsConformalAsPossibleVolumeParameterizer.h
)

set(SOURCES
  # Mappings
  Mapping.cc
  PiecewiseLinearMap.cc
  FundamentalMap.cc
  HarmonicFundamentalMap.cc
  BiharmonicFundamentalMap.cc
  # Surface mapping
  # ...
  # Volumetric mapping
  VolumeParameterizer.cc
  TetrahedralVolumeParameterizer.cc
  LinearTetrahedralVolumeParameterizer.cc
  FundamentalVolumeParameterizer.cc
  HarmonicTetrahedralVolumeParameterizer.cc
  HarmonicFundamentalVolumeParameterizer.cc
  AsConformalAsPossibleVolumeParameterizer.cc
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

# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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
  AsConformalAsPossibleVolumeParameterizer.h
  BiharmonicMap.h
  DiscreteMap.h
  FundamentalMap.h
  FundamentalVolumeParameterizer.h
  HarmonicFundamentalVolumeParameterizer.h
  HarmonicMap.h
  HarmonicTetrahedralVolumeParameterizer.h
  LinearTetrahedralVolumeParameterizer.h
  TetrahedralVolumeParameterizer.h
  VolumeParameterizer.h
  VolumetricMap.h
)

set(SOURCES
  AsConformalAsPossibleVolumeParameterizer.cc
  BiharmonicMap.cc
  DiscreteMap.cc
  FundamentalMap.cc
  FundamentalVolumeParameterizer.cc
  HarmonicFundamentalVolumeParameterizer.cc
  HarmonicMap.cc
  HarmonicTetrahedralVolumeParameterizer.cc
  LinearTetrahedralVolumeParameterizer.cc
  TetrahedralVolumeParameterizer.cc
  VolumeParameterizer.cc
  VolumetricMap.cc
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

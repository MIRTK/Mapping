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
# @brief Build configuration of MIRTK VolumetricMapping library.
##############################################################################

set(TARGET_NAME "Lib${PROJECT_NAME}")
set(OUTPUT_NAME "${PROJECT_PACKAGE_NAME}${PROJECT_NAME}")

set(HEADERS
  mirtkAsConformalAsPossibleVolumeParameterizer.h
  mirtkBiharmonicMap.h
  mirtkDiscreteMap.h
  mirtkFundamentalMap.h
  mirtkFundamentalVolumeParameterizer.h
  mirtkHarmonicFundamentalVolumeParameterizer.h
  mirtkHarmonicMap.h
  mirtkHarmonicTetrahedralVolumeParameterizer.h
  mirtkLinearTetrahedralVolumeParameterizer.h
  mirtkTetrahedralVolumeParameterizer.h
  mirtkVolumeParameterizer.h
  mirtkVolumetricMap.h
)

set(SOURCES
  mirtkAsConformalAsPossibleVolumeParameterizer.cc
  mirtkBiharmonicMap.cc
  mirtkDiscreteMap.cc
  mirtkFundamentalMap.cc
  mirtkFundamentalVolumeParameterizer.cc
  mirtkHarmonicFundamentalVolumeParameterizer.cc
  mirtkHarmonicMap.cc
  mirtkHarmonicTetrahedralVolumeParameterizer.cc
  mirtkLinearTetrahedralVolumeParameterizer.cc
  mirtkTetrahedralVolumeParameterizer.cc
  mirtkVolumeParameterizer.cc
  mirtkVolumetricMap.cc
)

foreach (hdr IN LISTS HEADERS)
  list(APPEND SOURCES "${PROJECT_INCLUDE_DIR}/${hdr}")
endforeach ()

basis_add_library(${TARGET_NAME} "${SOURCES}")

basis_set_target_properties(
  ${TARGET_NAME} PROPERTIES
    VERSION     "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
    SOVERSION   "${PROJECT_SOVERSION}"
    OUTPUT_NAME "${OUTPUT_NAME}"
)

basis_target_link_libraries(${TARGET_NAME}
  LibCommon
  LibNumerics
  LibImage
  LibPointSet
  ${VTK_LIBRARIES}
)

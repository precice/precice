# This file sets the XSDK defaults if specified
# Not necessary for precice are:
#   2. Install location
#   7. Determine precision
#      XSDK_PRECISION (default is double) and precice uses double
#   8. Determine index size
#      XSDK_INDEX_SIZE (default is 64) and precice uses size_t
#   9. Location of BLAS and LAPACK - not needed

# 1. XSDK Mode
# let the caller specify whether to use the xsdk defaults by default
set(USE_XSDK_DEFAULTS OFF CACHE BOOL "Use XSDK defaults?")
if(USE_XSDK_DEFAULTS)
  message(STATUS "XSDK defaults are enabled.")
  # 4. Debug and Release
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "The build type." FORCE)
  # 5. Shared and static libraries
  set(BUILD_SHARED_LIBS YES CACHE BOOL "Whether to build shared libraries." FORCE)
endif()

# 6. Build interface for additional language
if (DEFINED XSDK_ENABLE_CXX AND NOT XSDK_ENABLE_CXX)
  message(FATAL_ERROR "Using precice without CXX is pointless.")
endif()

if (DEFINED XSDK_ENABLE_PYTHON)
  message(STATUS "XSDK override: XSDK_ENABLE_PYTHON")
  set(PRECICE_PYTHONACTION "${XSDK_ENABLE_PYTHON}" CACHE BOOL "" FORCE)
endif()

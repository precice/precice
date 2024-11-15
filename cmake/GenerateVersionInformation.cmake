#
# This script generates the version information used in versions.hpp.in
#

# Adds a NAME=VALUE pair to the version information
macro(precice_vi_add_value NAME VALUE)
  set(${NAME} "${VALUE}")# CACHE INTERNAL "")
endmacro()

# Adds a boolean option to the version information
macro(precice_vi_add_option NAME)
  set(PVI_VAL "N")
  if(${NAME})
    set(PVI_VAL "Y")
  endif()
  precice_vi_add_value("${NAME}_VAL" ${PVI_VAL})
  unset(PVI_VAL)
endmacro()

precice_vi_add_option(PRECICE_FEATURE_MPI_COMMUNICATION)
precice_vi_add_option(PRECICE_FEATURE_PETSC_MAPPING)
precice_vi_add_option(PRECICE_FEATURE_GINKGO_MAPPING)
precice_vi_add_option(PRECICE_FEATURE_PYTHON_ACTIONS)
precice_vi_add_option(PRECICE_BINDINGS_C)
precice_vi_add_option(PRECICE_BINDINGS_FORTRAN)

# Add compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  precice_vi_add_value(PRECICE_VI_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  precice_vi_add_value(PRECICE_VI_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  precice_vi_add_value(PRECICE_VI_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
else()
  precice_vi_add_value(PRECICE_VI_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# Add linker flags
if(BUILD_SHARED_LIBS)
  precice_vi_add_value(PRECICE_VI_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
else()
  precice_vi_add_value(PRECICE_VI_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS}")
endif()

#
# This script generates the version information used in versions.hpp.in
#

# We use an internal cache to keep the code below a bit more sane
set(preCICE_VERSION_INFORMATION "" CACHE INTERNAL "")

# Adds a NAME=VALUE pair to the version information
function(precice_vi_add_value NAME VALUE)
  if(preCICE_VERSION_INFORMATION)
    # We have to use \; as separator to prevent the list-escaping of cmake
    # to remove the ; and escape the whitespaces
    set(preCICE_VERSION_INFORMATION "${preCICE_VERSION_INFORMATION}\;${NAME}=${VALUE}" CACHE INTERNAL "")
  else()
    set(preCICE_VERSION_INFORMATION "${NAME}=${VALUE}" CACHE INTERNAL "")
  endif()
endfunction()

# Adds a boolean option to the version information
function(precice_vi_add_option NAME)
  set(PVI_VAL "N")
  if(${NAME})
    set(PVI_VAL "Y")
  endif()
  precice_vi_add_value(${NAME} ${PVI_VAL})
endfunction()

precice_vi_add_option(PRECICE_MPICommunication)
precice_vi_add_option(PRECICE_PETScMapping)
precice_vi_add_option(PRECICE_PythonActions)
precice_vi_add_option(PRECICE_ENABLE_C)
precice_vi_add_option(PRECICE_ENABLE_FORTRAN)
precice_vi_add_value(CXX "${CMAKE_CXX_COMPILER_ID}")

# Add compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "")
  precice_vi_add_value(CXXFLAGS "${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  precice_vi_add_value(CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  precice_vi_add_value(CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  precice_vi_add_value(CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

# Add linker flags
if(BUILD_SHARED_LIBS)
  precice_vi_add_value(LDFAGS "${CMAKE_SHARED_LINKER_FLAGS}")
else()
  precice_vi_add_value(LDFAGS "${CMAKE_STATIC_LINKER_FLAGS}")
endif()

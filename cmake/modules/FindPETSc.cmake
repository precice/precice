#
# FindPETSc.cmake
#
# PETSc_NO_PATH - Runs in a clean environment
# PETSc_NO_ENV  - Do not use PETSc_DIR and PETSc_ARCH
#


#
# This function resolves a list of libraries based on a list of directories.
# Use this function to resolve the paths of a list of library names using a list of library directories.
# This is especially useful to handle the result of PkgConfig
#
function(_resolve_library_locations PREFIX LIBS DIRS)  
  unset(_libs)
  foreach(lib ${LIBS})
    find_library(${PREFIX}_LIB_${lib}
      NAMES ${lib}
      PATHS ${DIRS}
      NO_DEFAULT_PATH
      )
    list(APPEND _libs ${${PREFIX}_LIB_${lib}})
  endforeach(lib)
  set(${PREFIX}_LIBRARIES_RESOLVED "${_libs}" PARENT_SCOPE)
endfunction(_resolve_library_locations)


find_package(PkgConfig REQUIRED QUIET)

set(_env_path "")
if(NOT PETSc_NO_ENV AND DEFINED ENV{PETSC_DIR})
    if(DEFINED ENV{PETSC_ARCH})
      message(STATUS "Detected environment PETSC_DIR and PETSC_ARCH. Set PETSc_NO_ENV to disable this.")
      set(_env_path "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
    else()
      message(STATUS "Detected ENV PETSC_DIR without PETSC_ARCH. Set PETSc_NO_ENV to disable this.")
      set(_env_path "$ENV{PETSC_DIR}")
    endif()
endif()
list(APPEND CMAKE_PREFIX_PATH "${_env_path}")
unset(_env_path)

set(_pkg_args "")
if(PETSc_NO_PATH)
  set(_pkg_args "NO_CMAKE_ENVIRONMENT_PATH")
endif()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.5") 
  pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET GLOBAL ${_pkg_args} PETSc>=3.6)
  if(NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc ALIAS PkgConfig::PETSc)
  endif()
else()
  pkg_check_modules(PETSc REQUIRED ${args} PETSc>=3.6)
  if(NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc INTERFACE IMPORTED GLOBAL)
    unset(_libs)
    _resolve_library_locations(PETSc "${PETSc_LIBRARIES}" "${PETSc_LIBRARY_DIRS}")
    set_target_properties(PETSc::PETSc PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${PETSc_LIBRARIES_RESOLVED}"
      INTERFACE_COMPILE_OPTIONS "${PETSc_CFLAGS_OTHER}"
      ) 
  endif()
endif()

set(_versions "")
set(PETSc_VERSION_MAJOR "")
set(PETSc_VERSION_MINOR "")
string(REGEX MATCHALL "[0-9]+" _versions ${PETSc_VERSION})
list(GET _versions 0 PETSc_VERSION_MAJOR)
list(GET _versions 1 PETSc_VERSION_MINOR)

unset(_pkg_args)
unset(_versions)

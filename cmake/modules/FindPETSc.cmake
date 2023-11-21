# FindPETSc
# ---------
#
# Locates the PETSc library using pkg-config module PETSc
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following IMPORTED target:
#
#  PETSc::PETSc        - the PETSc library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
#  PETSc_FOUND          - if false, do not try to link to PETSc
#  PETSc_LIBRARIES      - a list of the full paths to all libraries
#  PETSc_INCLUDE_DIRS   - a list of all include directories
#  PETSc_VERSION        - the full version of PETSc MAJOR.MINOR.PATCH
#  PETSc_VERSION_MAJOR  - the MAJOR part of PETSc_VERSION
#  PETSc_VERSION_MINOR  - the MINOR part of PETSc_VERSION
#  PETSc_VERSION_PATCH  - the PATCH part of PETSc_VERSION
#
# Author: Frédéric Simonis @fsimonis

cmake_policy(VERSION 3.10)

# Generate a argument for cmake pkg-config call
if(PETSc_FIND_QUIETLY)
  find_package(PkgConfig QUIET)
else()
  find_package(PkgConfig)
endif()

if(PKG_CONFIG_FOUND)
  # Build the pkg-config version spec
  set(_pkg_version_spec "")
  if(DEFINED PETSc_FIND_VERSION)
    if(PETSc_FIND_VERSION_EXACT)
      set(_pkg_version_spec "=${PETSc_FIND_VERSION}")
    else()
      set(_pkg_version_spec ">=${PETSc_FIND_VERSION}")
    endif()
  endif()

  # Allow system flags
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_LIBS} 1)

  # Use pkg-config to find PETSc
  set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")

  if(PETSc_FIND_QUIETLY)
    pkg_check_modules(PETSc QUIET IMPORTED_TARGET GLOBAL "PETSc${_pkg_version_spec}")
  else()
    pkg_check_modules(PETSc IMPORTED_TARGET GLOBAL "PETSc${_pkg_version_spec}")
  endif()

  unset(_pkg_version_spec)

  # Extract version parts from the version information
  if(PETSc_VERSION)
    set(_petsc_versions "")
    string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PETSc_VERSION})
    list(GET _petsc_versions 0 _petsc_version_major)
    list(GET _petsc_versions 1 _petsc_version_minor)
    list(GET _petsc_versions 2 _petsc_version_patch)

    set(PETSc_VERSION ${PETSc_VERSION} CACHE STRING "Full version of PETSc")
    set(PETSc_VERSION_MAJOR ${_petsc_version_major} CACHE INTERNAL "Major version of PETSc")
    set(PETSc_VERSION_MINOR ${_petsc_version_minor} CACHE INTERNAL "Minor version of PETSc")
    set(PETSc_VERSION_PATCH ${_petsc_version_patch} CACHE INTERNAL "Patch version of PETSc")

    unset(_petsc_versions)
    unset(_petsc_version_major)
    unset(_petsc_version_minor)
    unset(_petsc_version_patch)
  endif()
endif()
unset(_petsc_quiet_arg)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  REQUIRED_VARS PETSc_FOUND PETSc_INCLUDE_DIRS PETSc_LIBRARIES
  VERSION_VAR PETSc_VERSION
  )

if(NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc ALIAS PkgConfig::PETSc)
endif()

mark_as_advanced(PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_VERSION_MAJOR PETSc_VERSION_MINOR PETSc_VERSION_PATCH VERSION_VAR PETSc_VERSION)

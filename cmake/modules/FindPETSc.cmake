# FindPETSc
# ---------
#
# Locates the PETSc library using pkg-config
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
# 
# This module defines the followwing IMPORTED target:
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
# Variables for locating PETSc
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Additional CMake variables for locating PETSc
#  PETSc_DIR            - the path to the root directory of PETSc 
#  PETSc_ARCH           - the PETSc architecture
#  PETSc_NO_ENV         - instructs the module not to use the environment variables 'PETSC_DIR' and 'PETSC_ENV' to find PETSc
#
# Environment Variables for locating PETSc
#  PETSC_DIR            - the path to the root directory of PETSc, part of the PETSc installation process
#  PETSC_ARCH           - the PETSc architecture, part of the PETSc installation process
#
# Author: Frédéric Simonis @fsimonis

cmake_policy(VERSION 3.10)


# Macro to print the search context used by pkg-config
macro(_petsc_print_pkg_env)
  if(NOT PETSc_FIND_QUIETLY)
    set(_env_mess "pkg-config will search the following paths:")
    if(DEFINED ENV{PKG_CONFIG_PATH})
      set(_env_mess "${_env_mess}\nPKG_CONFIG_PATH")
      string(REPLACE ":" ";" _env_pkg_list "$ENV{PKG_CONFIG_PATH}")
      foreach(p IN LISTS _env_pkg_list)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
      unset(_env_pkg_list)
    endif()
    if(CMAKE_PREFIX_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_PREFIX_PATH")
      foreach(p IN LISTS CMAKE_PREFIX_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    if(CMAKE_FRAMEWORK_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_FRAMEWORK_PATH")
      foreach(p IN LISTS CMAKE_FRAMEWORK_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    if(CMAKE_APPBUNDLE_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_APPBUNDLE_PATH")
      foreach(p IN LISTS CMAKE_APPBUNDLE_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    message(STATUS "${_env_mess}")
    unset(_env_mess)
  endif()
endmacro()


# Message macro which respects the QUIET arguemnt of the package
macro(_message)
  if(NOT PETSc_FIND_QUIETLY)
    message(${ARGV})
  endif()
endmacro()

set(_petsc_quiet_arg "")
if(PETSc_FIND_QUIETLY)
  set(_petsc_quiet_arg "QUIET")
endif()
find_package(PkgConfig ${_petsc_quiet_arg})

if(PKG_CONFIG_FOUND)
  # Detect additional pefix paths
  set(_petsc_detected_prefixes "")
  if(DEFINED PETSc_DIR)
    list(APPEND _petsc_detected_prefixes "${PETSc_DIR}")
    if(DEFINED PETSc_ARCH)
      list(APPEND _petsc_detected_prefixes "${PETSc_DIR}/${PETSc_ARCH}")
    endif()
  endif()

  if(DEFINED ENV{PETSC_DIR} AND NOT PETSc_NO_ENV)
    list(APPEND _petsc_detected_prefixes "$ENV{PETSC_DIR}")
    if(DEFINED ENV{PETSC_ARCH})
      list(APPEND _petsc_detected_prefixes "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
    endif()
  endif()
  list(REMOVE_DUPLICATES _petsc_detected_prefixes)

  set(_petsc_prefixes "")
  set(_petsc_skipped_prefixes "")
  _message(STATUS "Detecting additional PETSc prefixes")
  foreach(prefix IN LISTS _petsc_detected_prefixes )
    if(EXISTS "${prefix}/lib/pkgconfig")
      _message(STATUS "Detected ${prefix}")
      list(APPEND _petsc_prefixes "${prefix}")
    else()          
      list(APPEND _petsc_skipped_prefixes "${prefix}")
    endif()
  endforeach()
  if(_petsc_skipped_prefixes)
    _message(STATUS "Skipped the following invalid prefixes: ${_petsc_skipped_prefixes}")
  endif()
  unset(_petsc_detected_prefixes)

  # Remember the previous state of CMAKE_PREFIX_PATH
  set(_petsc_prefix_unset True)
  if(DEFINED CMAKE_PREFIX_PATH)
    set(_petsc_prefix_unset False)
    set(_petsc_prefix_old ${CMAKE_PREFIX_PATH})
  endif()

  list(APPEND CMAKE_PREFIX_PATH ${_petsc_prefixes})
  _petsc_print_pkg_env()

  # Build the pkg-config version spec
  set(_pkg_version_spec "")
  if(DEFINED PETSc_FIND_VERSION)
    if(PETSc_FIND_VERSION_EXACT)
      set(_pkg_version_spec "=${PETSc_FIND_VERSION}")
    else()
      set(_pkg_version_spec ">=${PETSc_FIND_VERSION}")
    endif()
  endif()

  # Use pkg-config to find PETSc
  set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")
  pkg_check_modules(PC_PETSc ${_petsc_quiet_arg} "PETSc${_pkg_version_spec}")
  unset(PKG_CONFIG_USE_CMAKE_PREFIX_PATH)
  unset(_pkg_version_spec)

  # Restore the previous state of CMAKE_PREFIX_PATH
  if(_petsc_prefix_unset)
    unset(CMAKE_PREFIX_PATH)
  else()
    set(CMAKE_PREFIX_PATH "${_petsc_prefix_old}")
  endif()

  # Set straight forward result variables
  set(PETSc_FOUND ${PC_PETSc_FOUND})
  set(PETSc_INCLUDE_DIRS ${PC_PETSc_INCLUDE_DIRS})

  # libm is always required
  set(_petsc_libs "m")
  set(_petsc_missing_libs "")

  # Find main PETSc libraries
  foreach(_next_lib IN LISTS PC_PETSc_LIBRARIES)
    find_library(_petsc_lib_${_next_lib} NAMES ${_next_lib} HINTS ${PC_PETSc_LIBRARY_DIRS})
    if(_petsc_lib_${_next_lib})
      list(APPEND _petsc_libs "${_petsc_lib_${_next_lib}}")
    else()
      list(APPEND _petsc_missing_libs "${_next_lib}")
    endif()
  endforeach()

  # Link against MPI if it is used.
  # This adds all required link directories.
  foreach(_next_lib IN LISTS PC_PETSc_STATIC_LIBRARIES)
    if(_next_lib STREQUAL "mpi")
      find_package(MPI ${_petsc_quiet_arg})
      if(MPI_FOUND)
        # Prefer to use the CXX dependencies if enabled otherwise use the C
        if(DEFINED CMAKE_CXX_COMPILER)
          list(APPEND _petsc_libs "MPI::MPI_CXX")
        else()
          enable_language(C)
          list(APPEND _petsc_libs "MPI::MPI_C")
        endif()
        break()
      else()
        list(APPEND _petsc_missing_libs "MPI")
      endif()
    endif()
  endforeach()

  # Check if everything was detected
  if(_petsc_missing_libs AND NOT PETSc_FIND_QUIETLY)
    message("The following libraries were not detected: ${_petsc_missing_libs}")
  elseif(NOT _petsc_missing_libs)
    # Set the visible variable. This will let the module to succeed.
    set(PETSc_LIBRARIES "${_petsc_libs}")
  endif()
  unset(_petsc_libs)
  unset(_petsc_missing_libs)

  # Extract version parts from the version information
  if(PC_PETSc_VERSION)
    set(_petsc_versions "")
    string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PC_PETSc_VERSION})
    list(GET _petsc_versions 0 _petsc_version_major)
    list(GET _petsc_versions 1 _petsc_version_minor)
    list(GET _petsc_versions 2 _petsc_version_patch)

    set(PETSc_VERSION ${PC_PETSc_VERSION} CACHE STRING "Full version of PETSc")
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
  add_library(PETSc::PETSc INTERFACE IMPORTED)
  set_target_properties(PETSc::PETSc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${PETSc_LIBRARIES}"
    )
endif()

mark_as_advanced(PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_VERSION_MAJOR PETSc_VERSION_MINOR PETSc_VERSION_PATCH VERSION_VAR PETSc_VERSION)

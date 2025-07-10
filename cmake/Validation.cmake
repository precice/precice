#
# preCICE Validation of Libraries
#
# This file contains macros which validate some libraries by trying to compile
# a minimalisitic executable.
# They report their activity as STATUS messages.
# A failing validation is reported as a FATAL_ERROR.
#
# Variables
#  PRECICE_ALWAYS_VALIDATE_LIBS - Always actively validates the libs.
#

# General function to validate a given library
# precice_validate_lib( <code-to-compile> NAME <name> [LINK_LIBRARIES <lib> ...] [COMPILE_DEFINITIONS <def> ... ]
#
# Caches on success to speed up development. Set PRECICE_ALWAYS_VALIDATE_LIBS to disable this behaviour.
#
function(precice_validate_lib ARG_CODE)
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 1 ARG
    ""
    "NAME"
    "LINK_LIBRARIES;COMPILE_DEFINITIONS")
  if(NOT ARG_NAME)
    message(FATAL_ERROR "Argument required: NAME")
  endif()
  if(NOT ARG_CODE)
    message(FATAL_ERROR "Argument required: CODE")
  endif()

  set(_filename validate${ARG_NAME}.cpp)
  set(_wdir ${CMAKE_CURRENT_BINARY_DIR}/validation)
  file(MAKE_DIRECTORY ${_wdir})
  set(_cache_success PRECICE_VALIDATE_${ARG_NAME}_SUCCESS)

  message(STATUS "Validating ${ARG_NAME}")
  file(WRITE ${_wdir}/${_filename} "${ARG_CODE}")

  if(${_cache_success} AND NOT PRECICE_ALWAYS_VALIDATE_LIBS)
    message(STATUS "Validating ${ARG_NAME} - success [cached]")
  else()
    unset(VAL_SUCCESS)
    unset(VAL_OUT)
    try_compile(VAL_SUCCESS
      ${_wdir}
      ${_wdir}/${_filename}
      COMPILE_DEFINITIONS ${ARG_COMPILE_DEFINITIONS}
      LINK_LIBRARIES ${ARG_LINK_LIBRARIES}
      OUTPUT_VARIABLE VAL_OUT
      CXX_STANDARD 11
      )
    if(NOT VAL_SUCCESS)
      message(FATAL_ERROR "Validating ${ARG_NAME} - failure\n\n${VAL_OUT}")
    else()
      message(STATUS "Validating ${ARG_NAME} - success")
      set(${_cache_success} Yes CACHE BOOL "Cached successful validation of ${ARG_NAME}." FORCE)
      mark_as_advanced(${_cache_success})
    endif()
  endif()
endfunction()

# Validation for LibPython
macro(precice_validate_libpython)
  precice_validate_lib(
    "#include <Python.h>\nint main() { return 0; } "
    NAME LibPython
    LINK_LIBRARIES Python3::Python
    )
endmacro()

# Validation for NumPy
macro(precice_validate_numpy)
  precice_validate_lib(
    "#include <Python.h>\n#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n#include <numpy/arrayobject.h>\n int main() { import_array1(0); return 0; } "
    NAME NumPy
    LINK_LIBRARIES Python3::NumPy Python3::Python
    )
endmacro()

# Validation for LibXML2
# We check for the header libxml/SAX2.h as we use it in preCICE
macro(precice_validate_libxml2)
  precice_validate_lib(
    "#include <libxml/SAX2.h>\nint main() { return 0; } "
  NAME LibXml2
  LINK_LIBRARIES LibXml2::LibXml2
  )
endmacro()

# Validation for Eigen
macro(precice_validate_eigen)
  precice_validate_lib(
    "#include <Eigen/Core>\nint main() { return 0; } "
    NAME Eigen
    LINK_LIBRARIES Eigen3::Eigen)
endmacro()

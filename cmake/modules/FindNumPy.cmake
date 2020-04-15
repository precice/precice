#.rst:
#
# Find the include directory for ``numpy/arrayobject.h``.
# This module exclusively uses the information from the numpy python module.
#
# This module sets the following variables:
#
# ``NumPy_FOUND``
#   True if NumPy was found.
# ``NumPy_INCLUDE_DIRS``
#   The include directories needed to use NumpPy.
# ``NumPy_VERSION``
#   The version of NumPy found.
#
# The module will also explicitly define one cache variable:
#
# ``NumPy_INCLUDE_DIR``
#

if(NOT NumPy_FOUND)
  set(_find_extra_args)
  if(NumPy_FIND_REQUIRED)
    list(APPEND _find_extra_args REQUIRED)
  endif()
  if(NumPy_FIND_QUIET)
    list(APPEND _find_extra_args QUIET)
  endif()
  find_package(PythonInterp ${_find_extra_args})
  find_package(PythonLibs ${_find_extra_args})

  if(PYTHON_EXECUTABLE)
    unset(_numpy_include_dir)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
      -c "import numpy; print(numpy.get_include())"
      OUTPUT_VARIABLE _numpy_include_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
      -c "import numpy; print(numpy.__version__)"
      OUTPUT_VARIABLE NumPy_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
    find_path(NumPy_INCLUDE_DIR
      numpy/arrayobject.h
      PATH "${_numpy_include_dir}"
      NO_DEFAULT_PATH
      )
    unset(_numpy_include_dir)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set NumPy_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( NumPy
    REQUIRED_VARS NumPy_INCLUDE_DIR
    VERSION_VAR NumPy_VERSION
    )

if(NumPy_FOUND)
    set(NumPy_INCLUDE_DIRS ${NumPy_INCLUDE_DIR})

    if(NOT TARGET NumPy::NumPy)
        add_library(NumPy::NumPy INTERFACE IMPORTED)
        set_property(TARGET NumPy::NumPy PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${NumPy_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(NumPy_INCLUDE_DIR)

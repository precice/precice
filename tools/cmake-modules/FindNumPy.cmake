#.rst:
#
# Find the include directory for ``numpy/arrayobject.h`` as well as other NumPy tools like ``conv-template`` and
# ``from-template``.
#
# This module sets the following variables:
#
# ``NumPy_FOUND``
#   True if NumPy was found.
# ``NumPy_INCLUDE_DIRS``
#   The include directories needed to use NumpPy.
# ``NumPy_VERSION``
#   The version of NumPy found.
# ``NumPy_CONV_TEMPLATE_EXECUTABLE``
#   Path to conv-template executable.
# ``NumPy_FROM_TEMPLATE_EXECUTABLE``
#   Path to from-template executable.
#
# The module will also explicitly define one cache variable:
#
# ``NumPy_INCLUDE_DIR``
#
# .. note::
#
#     To support NumPy < v0.15.0 where ``from-template`` and ``conv-template`` are not declared as entry points,
#     the module emulates the behavior of standalone executables by setting the corresponding variables with the
#     path the the python interpreter and the path to the associated script. For example:
#     ::
#
#         set(NumPy_CONV_TEMPLATE_EXECUTABLE /path/to/python /path/to/site-packages/numpy/distutils/conv_template.py CACHE STRING "Command executing conv-template program" FORCE)
#
#         set(NumPy_FROM_TEMPLATE_EXECUTABLE /path/to/python /path/to/site-packages/numpy/distutils/from_template.py CACHE STRING "Command executing from-template program" FORCE)
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

  find_program(NumPy_CONV_TEMPLATE_EXECUTABLE NAMES conv-template)
  find_program(NumPy_FROM_TEMPLATE_EXECUTABLE NAMES from-template)

  if(PYTHON_EXECUTABLE)
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

    # XXX This is required to support NumPy < v0.15.0. See note in module documentation above.
    if(NOT NumPy_CONV_TEMPLATE_EXECUTABLE)
      execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        -c "from numpy.distutils import conv_template; print(conv_template.__file__)"
        OUTPUT_VARIABLE _numpy_conv_template_file
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        )
      set(NumPy_CONV_TEMPLATE_EXECUTABLE "${PYTHON_EXECUTABLE}" "${_numpy_conv_template_file}" CACHE STRING "Command executing conv-template program" FORCE)
    endif()

    # XXX This is required to support NumPy < v0.15.0. See note in module documentation above.
    if(NOT NumPy_FROM_TEMPLATE_EXECUTABLE)
      execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        -c "from numpy.distutils import from_template; print(from_template.__file__)"
        OUTPUT_VARIABLE _numpy_from_template_file
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        )
      set(NumPy_FROM_TEMPLATE_EXECUTABLE "${PYTHON_EXECUTABLE}" "${_numpy_from_template_file}" CACHE STRING "Command executing from-template program" FORCE)
    endif()
  endif()
endif()

find_path(NumPy_INCLUDE_DIR
  numpy/arrayobject.h
  PATHS "${_numpy_include_dir}" "${PYTHON_INCLUDE_DIR}"
  PATH_SUFFIXES numpy/core/include
  )

set(NumPy_INCLUDE_DIRS ${NumPy_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set NumPy_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy
                                  REQUIRED_VARS
                                    NumPy_INCLUDE_DIR
                                    NumPy_CONV_TEMPLATE_EXECUTABLE
                                    NumPy_FROM_TEMPLATE_EXECUTABLE
                                  VERSION_VAR NumPy_VERSION
                                  )

mark_as_advanced(NumPy_INCLUDE_DIR)

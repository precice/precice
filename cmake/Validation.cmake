#
# preCICE Validation of Libraries
#
# This file contains macros which validate some libraries by trying to compile
# a minimalisitic exectutable.
# They report their activity as STATUS messages.
# A failing validation is reported as a FATAL_ERROR.
#

# Validation for LibPython
macro(precice_validate_libpython)
  message(STATUS "Validating LibPython")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/validateLibPython.cpp "#include <Python.h>\nint main() { return 0; } ")
  unset(VAL_SUCCESS)
  unset(VAL_OUT)
  try_compile(VAL_SUCCESS
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/validateLibPython.cpp
    COMPILE_DEFINITIONS "-I ${PYTHON_INCLUDE_DIRS}"
    LINK_LIBRARIES ${PYTHON_LIBRARIES}
    OUTPUT_VARIABLE VAL_OUT
    )
  if(NOT VAL_SUCCESS)
    message(FATAL_ERROR "Validating LibPython - failure\n\n${VAL_OUT}")
  else()
    message(STATUS "Validating LibPython - success")
  endif()
endmacro()

# Validation for NumPy
macro(precice_validate_numpy)
  message(STATUS "Validating NumPy")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/validateNumPy.cpp "#include <Python.h>\n#include <numpy/arrayobject.h>\nint main() { return 0; } ")
  unset(VAL_SUCCESS)
  unset(VAL_OUT)
  try_compile(VAL_SUCCESS
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/validateNumPy.cpp
    COMPILE_DEFINITIONS "-I ${PYTHON_INCLUDE_DIRS}"
    LINK_LIBRARIES NumPy::NumPy ${PYTHON_LIBRARIES}
    OUTPUT_VARIABLE VAL_OUT
    )
  if(NOT VAL_SUCCESS)
    message(FATAL_ERROR "Validating NumPy - failure\n\n${VAL_OUT}")
  else()
    message(STATUS "Validating NumPy - success")
  endif()

endmacro()

# Validation for LibXML2
# We check for the header libxml/SAX.h as we use it in preCICE
macro(precice_validate_libxml2)
  message(STATUS "Validating LibXML2")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/validateLibXml2.cpp "#include <libxml/SAX.h>\nint main() { return 0; } ")
  unset(VAL_SUCCESS)
  unset(VAL_OUT)
  try_compile(VAL_SUCCESS
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/validateLibXml2.cpp
    COMPILE_DEFINITIONS "-I ${LIBXML2_INCLUDE_DIR}"
    LINK_LIBRARIES ${LIBXML2_LIBRARIES}
    OUTPUT_VARIABLE VAL_OUT
    )
  if(NOT VAL_SUCCESS)
    message(FATAL_ERROR "Validating LibXML2 - failure\n\n${VAL_OUT}")
  else()
    message(STATUS "Validating LibXML2 - success")
  endif()
endmacro()

# Validation for Eigen
macro(precice_validate_eigen)
  message(STATUS "Validating Eigen")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/validateEigen.cpp "#include <Eigen/Core>\nint main() { return 0; } ")
  unset(VAL_SUCCESS)
  unset(VAL_OUT)
  try_compile(VAL_SUCCESS
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/validateEigen.cpp
    LINK_LIBRARIES Eigen3::Eigen
    OUTPUT_VARIABLE VAL_OUT
    )
  if(NOT VAL_SUCCESS)
    message(FATAL_ERROR "Validating Eigen - failure\n\n${VAL_OUT}")
  else()
    message(STATUS "Validating Eigen - success")
  endif()
endmacro()

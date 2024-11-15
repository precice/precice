# Adds the tests for the precice-profiling

# Single directory tests

set(_tests
  from-common-dir
  one-parallel-solver
  one-parallel-solver-missing-ranks
  one-parallel-solver-different-ids
  one-serial-damaged
  one-serial-solver
  one-solver-multiple-runs
  two-mixed-solvers
  two-parallel-solvers
  two-serial-solvers
  )

set(_precice_profiling "${CMAKE_CURRENT_LIST_DIR}/precice-profiling")
set(_test_script "${CMAKE_CURRENT_LIST_DIR}/tests/simple-test.cmake")
set(_test_root "${PROJECT_BINARY_DIR}/TestOutput")

# Check if polars is installed
set(_test_message_suffix "")
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import polars"
  RESULTS_VARIABLE PRECICE_PYTHON_POLARS_CHECK
  OUTPUT_QUIET ERROR_QUIET)
if(NOT "${PRECICE_PYTHON_POLARS_CHECK}" EQUAL "0")
  message(STATUS "The python package polars is not installed. Profiling tests will be limited. Install polars using your package manager or \"pip install polars\" to enable all profiling tests.")
  set(_test_message_suffix " - limited")
endif()

foreach(test IN LISTS _tests)
  set(_test_name precice.tools.profiling.${test})
  set(_test_location "${CMAKE_CURRENT_LIST_DIR}/tests/${test}")
  set(_test_wd "${_test_root}/tools.profiling.${test}")
  file(MAKE_DIRECTORY ${_test_wd})
  message(STATUS "Test ${_test_name}${_test_message_suffix}")
  add_test(NAME ${_test_name}
    COMMAND ${CMAKE_COMMAND} -DPROFILING_SCRIPT=${_precice_profiling} -DTEST_FOLDER=${_test_location} -DFOLDER_ARG=${_test_location} -P ${_test_script}
    WORKING_DIRECTORY ${_test_wd}
    )
endforeach(test)

# Multi directory tests

set(test from-solver-dirs)
set(_test_name precice.tools.profiling.${test})
set(_test_location "${CMAKE_CURRENT_LIST_DIR}/tests/${test}")
set(_test_wd "${PROJECT_BINARY_DIR}/TestOutput/tools.profiling.${test}")
file(MAKE_DIRECTORY ${_test_wd})
message(STATUS "Test ${_test_name}${_test_message_suffix}")
add_test(NAME ${_test_name}
  COMMAND ${CMAKE_COMMAND} -DPROFILING_SCRIPT=${_precice_profiling} -DTEST_FOLDER=${_test_location} -DFOLDER_ARG=${_test_location}/A\;${_test_location}/B -P ${_test_script}
  WORKING_DIRECTORY ${_test_wd}
  )

unset(_precice_profiling)
unset(_test_location)
unset(_test_message_suffix)
unset(_test_name)
unset(_test_root)
unset(_test_script)
unset(_test_wd)

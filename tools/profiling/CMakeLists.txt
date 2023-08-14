# Adds the tests for the precice-events

# Single directory tests

set(_tests
  from-common-dir
  one-parallel-solver
  one-parallel-solver-different-ids
  one-serial-damaged
  one-serial-solver
  one-solver-multiple-runs
  two-mixed-solvers
  two-parallel-solvers
  two-serial-solvers
  )

set(_precice_events "${CMAKE_CURRENT_LIST_DIR}/precice-events")
set(_test_script "${CMAKE_CURRENT_LIST_DIR}/tests/simple-test.cmake")
set(_test_root "${PROJECT_BINARY_DIR}/TestOutput")

# Install polars in a venv if available
if(NOT PRECICE_TEST_VENV)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import venv" RESULTS_VARIABLE _precice_has_venv ERROR_QUIET)
  if("${_precice_has_venv}" EQUAL 0)
    message("Installing polars in a python venv to allow for full profiling tests.")
    execute_process(COMMAND ${Python3_EXECUTABLE} -m venv --clear venv WORKING_DIRECTORY ${_test_root} OUTPUT_QUIET)
    execute_process(COMMAND ${_test_root}/venv/bin/pip install --upgrade pip polars OUTPUT_QUIET)

    execute_process(COMMAND ${_test_root}/venv/bin/python -c "import polars" RESULTS_VARIABLE _precice_has_polars ERROR_QUIET OUTPUT_QUIET)
    if ("${_precice_has_polars}" EQUAL 0)
      set(PRECICE_TEST_VENV "${_test_root}/venv" CACHE INTERNAL "Location of venv for preCICE tests")
    else()
      message(WARNING "Installing polars in a venv failed.")
    endif()
  else()
    message("Python venv not found. Please install python-venv or python3-venv for full profiling tests.")
  endif()
endif()

foreach(test IN LISTS _tests)
  set(_test_name precice.tools.events.${test})
  set(_test_location "${CMAKE_CURRENT_LIST_DIR}/tests/${test}")
  set(_test_wd "${_test_root}/tools.events.${test}")
  file(MAKE_DIRECTORY ${_test_wd})
  message(STATUS "Test ${_test_name}")
  add_test(NAME ${_test_name}
    COMMAND ${CMAKE_COMMAND} -DEVENTS_SCRIPT=${_precice_events} -DTEST_FOLDER=${_test_location} -DFOLDER_ARG=${_test_location} -P ${_test_script}
    WORKING_DIRECTORY ${_test_wd}
    )
  # Pass the venv to the script
  if (PRECICE_TEST_VENV)
    set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT "VIRTUAL_ENV=${PRECICE_TEST_VENV}")
  endif()
endforeach(test)

# Multi directory tests

set(test from-solver-dirs)
set(_test_name precice.tools.events.${test})
set(_test_location "${CMAKE_CURRENT_LIST_DIR}/tests/${test}")
set(_test_wd "${PROJECT_BINARY_DIR}/TestOutput/tools.events.${test}")
file(MAKE_DIRECTORY ${_test_wd})
message(STATUS "Test ${_test_name}")
add_test(NAME ${_test_name}
  COMMAND ${CMAKE_COMMAND} -DEVENTS_SCRIPT=${_precice_events} -DTEST_FOLDER=${_test_location} -DFOLDER_ARG=${_test_location}/A\;${_test_location}/B -P ${_test_script}
  WORKING_DIRECTORY ${_test_wd}
  )
if (PRECICE_TEST_VENV)
  set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT "VIRTUAL_ENV=${PRECICE_TEST_VENV}")
endif()


unset(_precice_events)

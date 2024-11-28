macro(start_file)
  set(script "")
endmacro()

function(add_command NAME TEST_NAME)
  set(args "")
  foreach(arg ${ARGN})
    if(arg MATCHES "[^-./:a-zA-Z0-9_]")
      string(APPEND args " [==[${arg}]==]")
    else()
      string(APPEND args " ${arg}")
    endif()
  endforeach()
  string(APPEND script "${NAME}(${TEST_NAME} ${args})\n")
  set(script "${script}" PARENT_SCOPE)
endfunction()

function(add_test NAME RANKS)
  if(NAME MATCHES "MPIPorts|MPISinglePorts" AND (PRECICE_MPI_VERSION MATCHES "Open MPI|Intel"))
    # Test is unsupported by this MPI implementation
    return()
  endif()

  string(REGEX REPLACE "[/<>]" "_" test_dir "${NAME}")
  set(test_dir "${TEST_DIR}/${test_dir}")

  if (NOT EXISTS test_dir)
    FILE(MAKE_DIRECTORY ${test_dir})
  endif()

  # Get Labels for general test grouping
  if(NAME MATCHES "^([^/]*)/.*$")
    string(REPLACE "Tests" "" _label "${CMAKE_MATCH_1}")
    list(APPEND labels "${_label}")
    unset(_label)
  endif()
  # Label common keywords
  foreach(kw Ginkgo Sockets MPIPorts MPISinglePorts MPI Remeshing)
    if(NAME MATCHES ${kw})
      list(APPEND labels ${kw})
    endif()
  endforeach()
  # Transform into digestable labels list (every test should be labeled)
  if(labels)
    string(TOLOWER "${labels}" labels)
    list(JOIN labels "\;" labels)
  endif()

  if(RANKS STREQUAL "1")
    add_command(add_test
      "[=[precice.${NAME}]=]"
      ${_TEST_EXECUTABLE} "--run_test=${NAME}" "--log_level=message"
    )
  else()
    # --map-by=:OVERSUBSCRIBE works for OpenMPI(4+5) and MPICH, but not for Intel which doesn't need a flag
    set(_oversubscribe_FLAG "--map-by;:OVERSUBSCRIBE")
    if(PRECICE_MPI_VERSION MATCHES "Intel")
      set(_oversubscribe_FLAG "")
    endif()

    add_command(add_test
      "[=[precice.${NAME}]=]"
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG}
      ${RANKS}
      ${_oversubscribe_FLAG}
      ${PRECICE_CTEST_MPI_FLAGS}
      ${MPIEXEC_PREFLAGS}
      ${_TEST_EXECUTABLE} "--run_test=${NAME}" "--log_level=message"
      ${MPIEXEC_POSTFLAGS}
    )
  endif()


  add_command(set_tests_properties
    "[=[precice.${NAME}]=]"
    PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2"
    TIMEOUT "30"
    WORKING_DIRECTORY "${test_dir}"
    LABELS "${labels}"
 )

 set(script "${script}" PARENT_SCOPE)
endfunction()

macro(flush_file)
  file(WRITE "${_TEST_FILE}" "${script}")
  set(flush_tests_MODE APPEND PARENT_SCOPE)
  unset(script)
endmacro()


function(discover_tests)
  cmake_parse_arguments(
    PARSE_ARGV
    0
    ""
    ""
    "TEST_EXECUTABLE;TEST_DIR;TEST_FILE"
    ""
  )

  execute_process(
    COMMAND "${_TEST_EXECUTABLE}" --list_units
    WORKING_DIRECTORY "${_TEST_DIR}"
    TIMEOUT 4
    OUTPUT_VARIABLE output
    RESULT_VARIABLE result
  )

  if(NOT ${result} EQUAL 0)
    string(REPLACE "\n" "\n    " output "${output}")
    message(FATAL_ERROR
      "Error running test executable.\n"
      "  Path: '${_TEST_EXECUTABLE}'\n"
      "  Working directory: '${_TEST_DIR}'\n"
      "  Result: ${result}\n"
      "  Output:\n"
      "    ${output}\n"
    )
  endif()

  # Turn into a CMake list
  string(REPLACE "\n" ";" output "${output}")

  start_file()
  foreach(test IN LISTS output)
    if (NOT test STREQUAL "")
      string(REGEX MATCH "(.+) ([0-9?])" _ "${test}")
      if (CMAKE_MATCH_2 STREQUAL "?")
        message(WARNING "Test ${CMAKE_MATCH_1} does not define a test setup!")
    else()
      add_test("${CMAKE_MATCH_1}" ${CMAKE_MATCH_2})
    endif()
    endif()
  endforeach()
  flush_file()
endfunction()


if(CMAKE_SCRIPT_MODE_FILE)
  discover_tests(
    TEST_EXECUTABLE "${TEST_EXECUTABLE}"
    TEST_DIR "${TEST_DIR}"
    TEST_FILE "${TEST_FILE}"
  )
endif()

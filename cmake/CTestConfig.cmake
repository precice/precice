#
# CTest
#

set(PRECICE_TEST_TIMEOUT_LONG 180 CACHE STRING "The timeout in seconds for longer tests.")
set(PRECICE_TEST_TIMEOUT_NORMAL 40 CACHE STRING "The timeout in seconds for normal tests.")
set(PRECICE_TEST_TIMEOUT_SHORT 20 CACHE STRING "The timeout in seconds for shorter tests.")

set(PRECICE_TEST_DIR "${preCICE_BINARY_DIR}/TestOutput")
mark_as_advanced(PRECICE_TEST_DIR)

# Solverdummies will be build in $PRECICE_SOLVERDUMMY_DIR/LANG
set(PRECICE_SOLVERDUMMY_DIR "${preCICE_BINARY_DIR}/Solverdummies")

include(CheckLanguage)
check_language(C)
check_language(Fortran)

# Detect the wrapper script that runs 2 solvers in parallel
set(PRECICE_TEST_WRAPPER_SCRIPT "")
if(UNIX)
  set(PRECICE_TEST_WRAPPER_SCRIPT "${preCICE_SOURCE_DIR}/cmake/runsolverdummies.sh")
else()
  message(STATUS "Running solverdummies on your system is not supported. We will ignore affected tests.")
endif()
mark_as_advanced(PRECICE_TEST_WRAPPER_SCRIPT)

function(add_precice_test_build_solverdummy PAT_LANG)
  # Turn language to lowercase
  string(TOLOWER ${PAT_LANG} PAT_LANG)

  # Locate the source directory
  set(PAT_SRC_DIR "${preCICE_SOURCE_DIR}/examples/solverdummies/${PAT_LANG}")
  if(NOT IS_DIRECTORY ${PAT_SRC_DIR})
    message(FATAL_ERROR "There is no solverdummy for language \"${PAT_LANG}\"")
  endif()

  # We always prefix our tests
  set(PAT_FULL_NAME "precice.solverdummy.build.${PAT_LANG}")

  # Make sure the required compiler is available
  if(PAT_LANG STREQUAL "fortran")
    if(NOT CMAKE_Fortran_COMPILER OR NOT PRECICE_BINDINGS_FORTRAN)
      message(STATUS "Test ${PAT_FULL_NAME} - skipped")
      return()
    endif()
  elseif(PAT_LANG STREQUAL "c")
    if(NOT CMAKE_C_COMPILER OR NOT PRECICE_BINDINGS_C)
      message(STATUS "Test ${PAT_FULL_NAME} - skipped")
      return()
    endif()
  endif()

  # Generate build directory
  set(PAT_BIN_DIR "${PRECICE_SOLVERDUMMY_DIR}/${PAT_LANG}")
  file(MAKE_DIRECTORY "${PAT_BIN_DIR}")

  # Add the actual test
  message(STATUS "Test ${PAT_FULL_NAME}")
  add_test(NAME ${PAT_FULL_NAME}
    COMMAND ${CMAKE_CTEST_COMMAND}
    --build-and-test ${PAT_SRC_DIR} ${PAT_BIN_DIR}
    --build-generator ${CMAKE_GENERATOR}
    --build-options -Dprecice_DIR=${preCICE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )

  # Setting properties
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    WORKING_DIRECTORY "${PAT_BIN_DIR}"
    FIXTURES_SETUP "${PAT_LANG}-solverdummy"
    LABELS "solverdummy"
    TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
    )
endfunction(add_precice_test_build_solverdummy)

function(add_precice_test_run_solverdummies PAT_LANG_A PAT_LANG_B)
  # Turn languages to lowercase
  string(TOLOWER ${PAT_LANG_A} PAT_LANG_A)
  string(TOLOWER ${PAT_LANG_B} PAT_LANG_B)

  # We always prefix our tests
  set(PAT_NAME "solverdummy.run.${PAT_LANG_A}-${PAT_LANG_B}")
  set(PAT_FULL_NAME "precice.${PAT_NAME}")

  if(NOT PRECICE_TEST_WRAPPER_SCRIPT)
      message(STATUS "Test ${PAT_FULL_NAME} - skipped")
      return()
  endif()

  # Make sure all required compilers are available
  foreach(_lang IN ITEMS ${PAT_LANG_A} ${PAT_LANG_B})
    if(_lang STREQUAL "fortran")
      if(NOT CMAKE_Fortran_COMPILER OR NOT PRECICE_BINDINGS_FORTRAN)
        message(STATUS "Test ${PAT_FULL_NAME} - skipped")
        return()
      endif()
    elseif(_lang STREQUAL "c")
      if(NOT CMAKE_C_COMPILER OR NOT PRECICE_BINDINGS_C)
        message(STATUS "Test ${PAT_FULL_NAME} - skipped")
        return()
      endif()
    endif()
  endforeach()

  # Locate the solverdummy config
  set(PAT_CONFIG "${preCICE_SOURCE_DIR}/examples/solverdummies/precice-config.xml")
  if(NOT EXISTS ${PAT_CONFIG})
    message(FATAL_ERROR "CMake was unable to locate the solverdummy config!")
  endif()

  # Locate binary dir of solverdummy A
  set(PAT_BIN_DIR_A "${PRECICE_SOLVERDUMMY_DIR}/${PAT_LANG_A}")
  if(NOT IS_DIRECTORY ${PAT_BIN_DIR_A})
    message(FATAL_ERROR "There is no configured solverdummy for language ${PAT_LANG_A}!")
  endif()

  # Locate binary dir of solverdummy B
  set(PAT_BIN_DIR_B "${PRECICE_SOLVERDUMMY_DIR}/${PAT_LANG_B}")
  if(NOT IS_DIRECTORY ${PAT_BIN_DIR_B})
    message(FATAL_ERROR "There is no configured solverdummy for language ${PAT_LANG_B}!")
  endif()

  # Generate run directory
  set(PAT_RUN_DIR "${PRECICE_TEST_DIR}/${PAT_NAME}")
  file(MAKE_DIRECTORY "${PAT_RUN_DIR}")

  # Add the actual test
  message(STATUS "Test ${PAT_FULL_NAME}")
  add_test(NAME ${PAT_FULL_NAME}
    COMMAND ${CMAKE_COMMAND}
    -D Python3_EXECUTABLE=${Python3_EXECUTABLE}
    -D WRAPPER=${PRECICE_TEST_WRAPPER_SCRIPT}
    -D CHECKER=${preCICE_SOURCE_DIR}/tools/profiling/validate-rank-files
    -D DUMMY_A=${PAT_BIN_DIR_A}/solverdummy
    -D DUMMY_B=${PAT_BIN_DIR_B}/solverdummy
    -D DUMMY_RUN_DIR=${PAT_RUN_DIR}
    -D DUMMY_CONFIG=${PAT_CONFIG}
    -P ${preCICE_SOURCE_DIR}/cmake/runsolverdummies.cmake
    )

  # Setting properties
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    WORKING_DIRECTORY "${PAT_RUN_DIR}"
    FIXTURES_REQUIRED "${PAT_LANG_A}-solverdummy;${PAT_LANG_B}-solverdummy"
    LABELS "solverdummy"
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    )
endfunction(add_precice_test_run_solverdummies)


enable_testing()

# Autodiscovery of CTest
if(NOT PRECICE_FEATURE_MPI_COMMUNICATION)
  message(STATUS "Unit and integrationtests require MPI to be enabled.")
else()

  # This file is automatically loaded by CTEST and needs to exist to prevent strange errors
  set(ctest_tests_file "${preCICE_BINARY_DIR}/ctest_tests.cmake")
  if(NOT EXISTS "${ctest_tests_file}")
    file(WRITE "${ctest_tests_file}" "")
  endif()
  set_property(DIRECTORY
    APPEND PROPERTY TEST_INCLUDE_FILES "${ctest_tests_file}"
  )

  # Custom command to generate the tests list using the testprecice binary
  add_custom_command(
    OUTPUT "${ctest_tests_file}"
    COMMAND "${CMAKE_COMMAND}"
    -D "TEST_EXECUTABLE=$<TARGET_FILE:testprecice>"
    -D "TEST_FILE=${ctest_tests_file}"
    -D "TEST_DIR=${PRECICE_TEST_DIR}"
    -D "PRECICE_MPI_VERSION=${PRECICE_MPI_VERSION}"
    -D "MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
    -D "MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
    -D "PRECICE_CTEST_MPI_FLAGS=${PRECICE_CTEST_MPI_FLAGS}"
    -D "MPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
    -D "MPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}"
    -P "${preCICE_SOURCE_DIR}/cmake/discover_tests.cmake"
    DEPENDS testprecice
    COMMENT "Generating list of tests"
    VERBATIM)

  # Custom target that forces the test list to be updated
  add_custom_target(precice-test-list ALL
    DEPENDS "${ctest_tests_file}"
  )
endif()

# Add solverdummy tests
add_precice_test_build_solverdummy(cpp)
add_precice_test_build_solverdummy(c)
add_precice_test_build_solverdummy(fortran)

add_precice_test_run_solverdummies(cpp cpp)
add_precice_test_run_solverdummies(c c)
add_precice_test_run_solverdummies(fortran fortran)

add_precice_test_run_solverdummies(cpp c)
add_precice_test_run_solverdummies(cpp fortran)
add_precice_test_run_solverdummies(c fortran)

# Add tests for precice-tools

if(PRECICE_BUILD_TOOLS)
  function(add_tools_test)
    cmake_parse_arguments(PARSE_ARGV 0 PAT "WILL_FAIL" "NAME;MATCH" "COMMAND")
    set(PAT_FULL_NAME "precice.tools.${PAT_NAME}")
    message(STATUS "Test ${PAT_FULL_NAME}")
    add_test(
      NAME ${PAT_FULL_NAME}
      COMMAND ${PAT_COMMAND}
      )
    set_tests_properties(${PAT_FULL_NAME} PROPERTIES TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT} LABELS "tools;bin")
    if(PAT_WILL_FAIL)
      set_tests_properties(${PAT_FULL_NAME} PROPERTIES WILL_FAIL YES)
    endif()
    if(PAT_MATCH)
      set_tests_properties(${PAT_FULL_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "${PAT_MATCH}")
    endif()
  endfunction()

  # precice-tools

  add_tools_test(
    NAME legacy.noarg
    COMMAND precice-tools
    WILL_FAIL)

  add_tools_test(
    NAME legacy.invalidcmd
    COMMAND precice-tools invalidcommand
    WILL_FAIL)

  add_tools_test(
    NAME legacy.version
    COMMAND precice-tools version
    MATCH "${CMAKE_PROJECT_VERSION}"
    )

  add_tools_test(
    NAME legacy.versionopt
    COMMAND precice-tools --version
    MATCH "${CMAKE_PROJECT_VERSION}"
    )

  add_tools_test(
    NAME legacy.markdown
    COMMAND precice-tools md
    MATCH "# precice-configuration"
    )

  add_tools_test(
    NAME legacy.xml
    COMMAND precice-tools xml
    MATCH "<!-- TAG precice-configuration"
    )

  add_tools_test(
    NAME legacy.dtd
    COMMAND precice-tools dtd
    MATCH "<!ELEMENT precice-configuration"
    )

  add_tools_test(
    NAME legacy.check.file
    COMMAND precice-tools check ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml
    )

  add_tools_test(
    NAME legacy.check.file+name
    COMMAND precice-tools check ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml SolverTwo
    )

  add_tools_test(
    NAME legacy.check.file+name+size
    COMMAND precice-tools check ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml SolverTwo 2
    )

  # precice-version

  add_tools_test(
    NAME version
    COMMAND precice-version
    MATCH "${CMAKE_PROJECT_VERSION}"
    )

  # precice-config-doc

  add_tools_test(
    NAME config-doc.noarg
    COMMAND precice-config-doc
    WILL_FAIL)

  add_tools_test(
    NAME config-doc.markdown
    COMMAND precice-config-doc md
    MATCH "# precice-configuration"
    )

  add_tools_test(
    NAME config-doc.xml
    COMMAND precice-config-doc xml
    MATCH "<!-- TAG precice-configuration"
    )

  add_tools_test(
    NAME config-doc.dtd
    COMMAND precice-config-doc dtd
    MATCH "<!ELEMENT precice-configuration"
    )

  # precice-config-validate

  add_tools_test(
    NAME config-validate.noarg
    COMMAND precice-config-validate
    WILL_FAIL)

  add_tools_test(
    NAME config-validate.missingfile
    COMMAND precice-config-validate ${PROJECT_SOURCE_DIR}/definitely/missing/file.xml
    MATCH ERROR:
    )

  add_tools_test(
    NAME config-validate.file
    COMMAND precice-config-validate ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml
    )

  add_tools_test(
    NAME config-validate.file+name
    COMMAND precice-config-validate ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml SolverTwo
    )

  add_tools_test(
    NAME config-validate.file+name+size
    COMMAND precice-config-validate ${PROJECT_SOURCE_DIR}/src/precice/tests/config-checker.xml SolverTwo 2
    )

  # Simple configuration tests

  function(precice_test_config_valid path)
    set(name "precice.config.${path}.valid")
    set(solver "")
    set(ranks "")

    if (ARGC GREATER 1)
      set(solver "${ARGV1}")
      set(name "${name}:${solver}")
    endif()

    if (ARGC GREATER 2)
      set(ranks "${ARGV2}")
      set(name "${name}@${ranks}")
    endif()
    add_test(NAME ${name}
      COMMAND precice-config-validate "${PROJECT_SOURCE_DIR}/tests/config/${path}" ${name} ${ranks}
      )
    set_tests_properties(${name}
      PROPERTIES
      TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
      LABELS "tools;bin;config;valid")
  endfunction()

  function(precice_test_config_invalid path expression)
    set(name "precice.config.${path}.invalid")
    set(solver "")
    set(ranks "")

    if (ARGC GREATER 2)
      set(solver "${ARGV2}")
      set(name "${name}:${solver}")
    endif()

    if (ARGC GREATER 3)
      set(ranks "${ARGV3}")
      set(name "${name}@${ranks}")
    endif()
    add_test(NAME ${name}
      COMMAND precice-config-validate "${PROJECT_SOURCE_DIR}/tests/config/${path}" ${name} ${ranks}
      )
    set_tests_properties(${name}
      PROPERTIES
      TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
      PASS_REGULAR_EXPRESSION "${expression}"
      LABELS "tools;bin;config;invalid")
  endfunction()

  include(${PROJECT_SOURCE_DIR}/tests/config/tests.cmake)
endif()

# Add a separate target to test only the base
add_custom_target(
  test_base
  COMMAND ctest -V
  DEPENDS testprecice
  WORKING_DIRECTORY ${preCICE_BINARY_DIR}
  USES_TERMINAL
  )

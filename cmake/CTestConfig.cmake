#
# CTest
#

set(PRECICE_TEST_TIMEOUT_LONG 60 CACHE STRING "The timeout in seconds for longer tests.")
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


function(add_precice_test)
  cmake_parse_arguments(PARSE_ARGV 0 PAT "PETSC;MPIPORTS" "NAME;ARGUMENTS;TIMEOUT;LABELS" "")
  # Check arguments
  if(NOT PAT_NAME)
    message(FATAL_ERROR "Argument NAME not passed")
  endif()

  # We always prefix our tests
  set(PAT_FULL_NAME "precice.${PAT_NAME}")

  # Are direct dependencies fullfilled?
  if( (NOT PRECICE_MPICommunication) OR (PAT_PETSC AND NOT PRECICE_PETScMapping) )
    message(STATUS "Test ${PAT_FULL_NAME} - skipped")
    return()
  endif()

  if(PAT_MPIPORTS AND PRECICE_MPI_OPENMPI)
    message(STATUS "Test ${PAT_FULL_NAME} - skipped (OpenMPI)")
    return()
  endif()

  # Assemble the command
  message(STATUS "Test ${PAT_FULL_NAME}")
  add_test(NAME ${PAT_FULL_NAME}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${PRECICE_CTEST_MPI_FLAGS} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:testprecice> ${MPIEXEC_POSTFLAGS} ${PAT_ARGUMENTS}
    )
  # Generate working directory
  set(PAT_WDIR "${PRECICE_TEST_DIR}/${PAT_NAME}")
  file(MAKE_DIRECTORY "${PAT_WDIR}")
  # Setting properties
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    WORKING_DIRECTORY "${PAT_WDIR}"
    ENVIRONMENT "PRECICE_ROOT=${preCICE_SOURCE_DIR};OMPI_MCA_rmaps_base_oversubscribe=1"
    )
  if(PAT_TIMEOUT)
    set_tests_properties(${PAT_FULL_NAME} PROPERTIES TIMEOUT ${PAT_TIMEOUT} )
  endif()
  set(_labels ${PAT_LABELS})
  set_tests_properties(${PAT_FULL_NAME} PROPERTIES LABELS "${_labels}")
endfunction(add_precice_test)

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
    if(NOT CMAKE_Fortran_COMPILER OR NOT PRECICE_ENABLE_FORTRAN)
      message(STATUS "Test ${PAT_FULL_NAME} - skipped")
      return()
    endif()
  elseif(PAT_LANG STREQUAL "c")
    if(NOT CMAKE_C_COMPILER OR NOT PRECICE_ENABLE_C)
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
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    WORKING_DIRECTORY "${PAT_BIN_DIR}"
    FIXTURES_SETUP "${PAT_LANG}-solverdummy"
    LABELS "Solverdummy"
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
      if(NOT CMAKE_Fortran_COMPILER OR NOT PRECICE_ENABLE_FORTRAN)
        message(STATUS "Test ${PAT_FULL_NAME} - skipped")
        return()
      endif()
    elseif(_lang STREQUAL "c")
      if(NOT CMAKE_C_COMPILER OR NOT PRECICE_ENABLE_C)
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
    -D WRAPPER=${PRECICE_TEST_WRAPPER_SCRIPT}
    -D DUMMY_A=${PAT_BIN_DIR_A}/solverdummy
    -D DUMMY_B=${PAT_BIN_DIR_B}/solverdummy
    -D DUMMY_RUN_DIR=${PAT_RUN_DIR}
    -D DUMMY_CONFIG=${PAT_CONFIG}
    -P ${preCICE_SOURCE_DIR}/cmake/runsolverdummies.cmake
    )

  # Setting properties
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    WORKING_DIRECTORY "${PAT_RUN_DIR}"
    FIXTURES_REQUIRED "${PAT_LANG_A}-solverdummy"
    FIXTURES_REQUIRED "${PAT_LANG_B}-solverdummy"
    LABELS "Solverdummy"
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    )
endfunction(add_precice_test_run_solverdummies)


enable_testing()

if(NOT PRECICE_MPICommunication)
  message("Tests require MPICommunication to be enabled.")
endif()

add_precice_test(
  NAME acceleration
  ARGUMENTS "--run_test=AccelerationTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME action
  ARGUMENTS "--run_test=ActionTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME com
  ARGUMENTS "--run_test=CommunicationTests:\!CommunicationTests/MPIPorts:\!CommunicationTests/MPISinglePorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME com.mpiports
  ARGUMENTS "--run_test=CommunicationTests/MPIPorts:CommunicationTests/MPISinglePorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  LABELS "mpiports"
  MPIPORTS
  )
add_precice_test(
  NAME cplscheme
  ARGUMENTS "--run_test=CplSchemeTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME io
  ARGUMENTS "--run_test=IOTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME m2n
  ARGUMENTS "--run_test=M2NTests:\!M2NTests/MPIPorts:\!M2NTets/MPISinglePorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME m2n.mpiports
  ARGUMENTS "--run_test=M2NTests/MPIPorts:M2NTests/MPISinglePorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  LABELS "mpiports"
  MPIPORTS
  )
add_precice_test(
  NAME mapping
  ARGUMENTS "--run_test=MappingTests:\!MappingTests/PetRadialBasisFunctionMapping"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME mapping.petrbf
  ARGUMENTS "--run_test=MappingTests/PetRadialBasisFunctionMapping"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  LABELS petsc
  PETSC
  )
add_precice_test(
  NAME math
  ARGUMENTS "--run_test=MathTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME mesh
  ARGUMENTS "--run_test=MeshTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME partition
  ARGUMENTS "--run_test=PartitionTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME interface
  ARGUMENTS "--run_test=PreciceTests:\!PreciceTests/Serial:\!PreciceTests/Parallel"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME serial
  ARGUMENTS "--run_test=PreciceTests/Serial"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
  )
add_precice_test(
  NAME parallel
  ARGUMENTS "--run_test=PreciceTests/Parallel"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
  )
add_precice_test(
  NAME query
  ARGUMENTS "--run_test=QueryTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME testing
  ARGUMENTS "--run_test=TestingTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME utils
  ARGUMENTS "--run_test=UtilsTests"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME xml
  ARGUMENTS "--run_test=XML"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )

add_precice_test_build_solverdummy(cpp)
add_precice_test_build_solverdummy(c)
add_precice_test_build_solverdummy(fortran)

add_precice_test_run_solverdummies(cpp cpp)
add_precice_test_run_solverdummies(c c)
add_precice_test_run_solverdummies(fortran fortran)

add_precice_test_run_solverdummies(cpp c)
add_precice_test_run_solverdummies(cpp fortran)
add_precice_test_run_solverdummies(c fortran)

# Add a separate target to test only the base
add_custom_target(
  test_base
  COMMAND ctest -V
  DEPENDS testprecice
  WORKING_DIRECTORY ${preCICE_BINARY_DIR}
  USES_TERMINAL
  )

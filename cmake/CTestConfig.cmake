#
# CTest
#

set(PRECICE_TEST_TIMEOUT_LONG 180 CACHE STRING "The timeout in seconds for longer tests.")
set(PRECICE_TEST_TIMEOUT_SHORT 20 CACHE STRING "The timeout in seconds for shorter tests.")

set(PRECICE_TEST_DIR "${preCICE_BINARY_DIR}/TestOutput")
mark_as_advanced(PRECICE_TEST_DIR)

function(add_precice_test)
  cmake_parse_arguments(PARSE_ARGV 0 PAT "MPI;CANFAIL" "NAME;ARGUMENTS;TIMEOUT" "")
  if(NOT PAT_NAME)
    message(FATAL_ERROR "Argument NAME not passed")
  endif()
  # We always prefix our tests
  set(PAT_FULL_NAME "precice.${PAT_NAME}")
  message(STATUS "Adding Test ${PAT_FULL_NAME}")
  # Generate working directory
  set(PAT_WDIR "${PRECICE_TEST_DIR}/${PAT_NAME}")
  file(MAKE_DIRECTORY "${PAT_WDIR}")
  # Assemble the command

  if(PAT_MPI)
    add_test(NAME ${PAT_FULL_NAME}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${PRECICE_CTEST_MPI_FLAGS} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS} --log_level=all ${MPIEXEC_POSTFLAGS}
      )
  else()
    add_test(NAME ${PAT_FULL_NAME}
      COMMAND $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS} --log_level=all
      )
  endif()

  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    WORKING_DIRECTORY "${PAT_WDIR}"
    ENVIRONMENT PRECICE_ROOT=${preCICE_SOURCE_DIR}
    )
  if(PAT_TIMEOUT)
    set_tests_properties(${PAT_FULL_NAME} PROPERTIES TIMEOUT ${PAT_TIMEOUT} )
  endif()
  if(PAT_CANFAIL)
    set_tests_properties(${PAT_FULL_NAME} PROPERTIES LABELS "canfail")
  endif()
endfunction(add_precice_test)

enable_testing()

if(MPI AND MPIEXEC_EXECUTABLE)
  # Configured with MPI
  # Register the uncritical base-set
  add_precice_test(
    NAME Base
    ARGUMENTS "--run_test=\!@MPI_Ports:\!MappingTests/PetRadialBasisFunctionMapping/Parallel/\*"
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    MPI
    )
  add_precice_test(
    NAME MPI_Ports
    ARGUMENTS "--run_test=@MPI_Ports"
    TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}

    MPI
    CANFAIL
    )
  if(PETSC)
    add_precice_test(
      NAME PetRBFParallel
      ARGUMENTS "--run_test=MappingTests/PetRadialBasisFunctionMapping/Parallel/\*"
      TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
      MPI
      CANFAIL
      )
  endif()
  add_precice_test(
    NAME NoMPI
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    )
else()
  # Configured without MPI
  add_precice_test(
    NAME Base
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    )
endif()


# Add a separate target to test only the base
add_custom_target(
  test_base
  COMMAND ctest -V -R precice.Base
  DEPENDS testprecice
  WORKING_DIRECTORY ${preCICE_BINARY_DIR}
  USES_TERMINAL
  )

#
# CTest
#

set(PRECICE_TEST_DIR "${preCICE_BINARY_DIR}/TestOutput")
mark_as_advanced(PRECICE_TEST_DIR)

function(add_precice_test)
  cmake_parse_arguments(PARSE_ARGV 0 PAT "MPI;CANFAIL" "NAME;ARGUMENTS" "")
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
      COMMAND ${MPIEXEC_EXECUTABLE} -np 4 $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS}
      )
  else()
    add_test(NAME ${PAT_FULL_NAME}
      COMMAND $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS}
      )
  endif()
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    TIMEOUT 20 # Set the timeout to 60 seconds on this test
    WORKING_DIRECTORY "${PAT_WDIR}"
    )
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
    MPI
    )
  add_precice_test(
    NAME MPI_Ports
    ARGUMENTS "--run_test=@MPI_Ports"
    MPI
    CANFAIL
    )
  add_precice_test(
    NAME PetRBFParallel
    ARGUMENTS "--run_test=MappingTests/PetRadialBasisFunctionMapping/Parallel/\*"
    MPI
    CANFAIL
    )
  add_precice_test(
    NAME Full
    MPI
    CANFAIL
    )
  add_precice_test(
    NAME NoMPI
    )
else()
  # Configured without MPI
  add_precice_test(
    NAME NoMPI
    )
endif()

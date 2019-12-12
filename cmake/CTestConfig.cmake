#
# CTest
#

set(PRECICE_TEST_TIMEOUT_LONG 180 CACHE STRING "The timeout in seconds for longer tests.")
set(PRECICE_TEST_TIMEOUT_SHORT 20 CACHE STRING "The timeout in seconds for shorter tests.")

set(PRECICE_TEST_DIR "${preCICE_BINARY_DIR}/TestOutput")
mark_as_advanced(PRECICE_TEST_DIR)

function(add_precice_test)
  cmake_parse_arguments(PARSE_ARGV 0 PAT "NOMPI;PETSC;MPI;CANFAIL" "NAME;ARGUMENTS;TIMEOUT;LABELS" "")
  # Check arguments
  if(NOT PAT_NAME)
    message(FATAL_ERROR "Argument NAME not passed")
  endif()
  if(PAT_MPI AND PAT_NOMPI)
    message(FATAL_ERROR "You cannot specify MPI and NOMPI simultaneously.")
  endif()
  if(PAT_NOMPI AND PAT_PETSC)
    message(FATAL_ERROR "You cannot specify NOMPI and PETSC simultaneously.")
  endif()

  # We always prefix our tests
  set(PAT_FULL_NAME "precice.${PAT_NAME}")

  # Are direct dependencies fullfilled?
  if( (PAT_MPI AND NOT MPI) OR (PAT_PETSC AND NOT PETSC) )
    message(STATUS "Test ${PAT_FULL_NAME} - skipped")
    return()
  endif()

  # Assemble the command
  if(PAT_PETSC OR (NOT PAT_NOMPI AND MPI))
    # Parallel tests, dispatched by MPI
    message(STATUS "Test ${PAT_FULL_NAME} - parallel")
    add_test(NAME ${PAT_FULL_NAME}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${PRECICE_CTEST_MPI_FLAGS} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS} ${MPIEXEC_POSTFLAGS}
      )
  elseif(NOT PAT_MPI)
    # Serial tests, called directly
    message(STATUS "Test ${PAT_FULL_NAME} - serial")
    add_test(NAME ${PAT_FULL_NAME}
      COMMAND $<TARGET_FILE:testprecice> ${PAT_ARGUMENTS}
      )
  else()
    message(STATUS "Test ${PAT_FULL_NAME} - skipped")
    return()
  endif()
  # Generate working directory
  set(PAT_WDIR "${PRECICE_TEST_DIR}/${PAT_NAME}")
  file(MAKE_DIRECTORY "${PAT_WDIR}")
  # Setting properties
  set_tests_properties(${PAT_FULL_NAME}
    PROPERTIES
    RUN_SERIAL TRUE # Do not run this test in parallel with others
    WORKING_DIRECTORY "${PAT_WDIR}"
    ENVIRONMENT PRECICE_ROOT=${preCICE_SOURCE_DIR}
    )
  if(PAT_TIMEOUT)
    set_tests_properties(${PAT_FULL_NAME} PROPERTIES TIMEOUT ${PAT_TIMEOUT} )
  endif()
  set(_labels ${PAT_LABELS})
  if(PAT_CANFAIL)
    list(APPEND _labels "canfail")
  endif()
  set_tests_properties(${PAT_FULL_NAME} PROPERTIES LABELS "${_labels}")
endfunction(add_precice_test)

enable_testing()

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
  ARGUMENTS "--run_test=CommunicationTests:\!CommunicationTests/MPIPorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  MPI
  )
add_precice_test(
  NAME com.mpiports
  ARGUMENTS "--run_test=CommunicationTests/MPIPorts"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  LABELS "mpiports;canfail"
  MPI
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
  ARGUMENTS "--run_test=M2NTests:\!M2NTests/MPIPortsCommunication"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  MPI
  )
add_precice_test(
  NAME m2n.mpiports
  ARGUMENTS "--run_test=M2NTests/MPIPortsCommunication"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  LABELS "mpiports;canfail"
  MPI
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
  MPI
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
  MPI
  )
add_precice_test(
  NAME interface
  ARGUMENTS "--run_test=PreciceTests:\!PreciceTests/Serial:\!PreciceTests/Parallel"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  )
add_precice_test(
  NAME serial
  ARGUMENTS "--run_test=PreciceTests/Serial"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  MPI
  )
add_precice_test(
  NAME parallel
  ARGUMENTS "--run_test=PreciceTests/Parallel"
  TIMEOUT ${PRECICE_TEST_TIMEOUT_SHORT}
  MPI
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

if(MPI)
  add_precice_test(
    NAME nompi
    TIMEOUT ${PRECICE_TEST_TIMEOUT_LONG}
    NOMPI
    )
endif()

# Add a separate target to test only the base
add_custom_target(
  test_base
  COMMAND ctest -V -LE canfail
  DEPENDS testprecice
  WORKING_DIRECTORY ${preCICE_BINARY_DIR}
  USES_TERMINAL
  )

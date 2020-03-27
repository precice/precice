if(NOT EXISTS ${DUMMY_A})
  message(FATAL_ERROR "CMake was unable to locate solverdummy A at ${DUMMY_A}!")
endif()

if(NOT EXISTS ${DUMMY_B})
  message(FATAL_ERROR "CMake was unable to locate solverdummy A at ${DUMMY_B}!")
endif()

if(NOT EXISTS ${DUMMY_RUN_DIR})
  message(FATAL_ERROR "CMake was unable to locate the working directory at ${DUMMY_RUN_DIR}!")
endif()

if(NOT EXISTS ${DUMMY_CONFIG})
  message(FATAL_ERROR "CMake was unable to locate the solverdummy config at ${DUMMY_CONFIG}!")
endif()

# Remove precice-run
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory "${DUMMY_RUN_DIR}/precice-run")

# Execute both solvers in parallel
execute_process(
  COMMAND ${DUMMY_A} ${DUMMY_CONFIG} SolverOne MeshOne
  COMMAND ${DUMMY_B} ${DUMMY_CONFIG} SolverTwo MeshTwo
  WORKING_DIRECTORY ${DUMMY_RUN_DIR}
  RESULTS_VARIABLE DUMMY_RESULTS
)

# Check the return codes of the solvers
foreach(RE IN LISTS DUMMY_RESULTS)
  if(NOT RE EQUAL 0)
    message(FATAL_ERROR "An error occured running the solverdummies!")
  endif()
endforeach(RE)

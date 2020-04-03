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

# Try at most 10 times
foreach(_current_try RANGE 1 10)

  # Execute both solvers in parallel
  execute_process(
    COMMAND ${DUMMY_A} ${DUMMY_CONFIG} SolverOne MeshOne
    COMMAND ${DUMMY_B} ${DUMMY_CONFIG} SolverTwo MeshTwo
    WORKING_DIRECTORY ${DUMMY_RUN_DIR}
    RESULTS_VARIABLE DUMMY_RESULTS
  )

  # Check the return codes/statuses of the solvers
  set(_rerun_test False)
  foreach(RE IN LISTS DUMMY_RESULTS)
    if((RE STREQUAL "SIGPIPE"))
      # Rerun the test in case we encounter a SIGPIPE
      message(STATUS "The test stopped due to SIGPIPE. Rerunning ...")
      set(_rerun_test True)
    elseif(NOT (RE EQUAL 0))
      # Fail in case we encounter another error code/condition other than 0
      message(FATAL_ERROR "An error occured running the solverdummies! Return codes : \"${DUMMY_RESULTS}\"")
    endif()
  endforeach(RE)

  # Only rerun if we need to.
  if(_rerun_test)
    if(_current_try EQUAL 10)
      message(FATAL_ERROR "The test failed due to SIGPIPE 10 times in a row.")
    endif()
  else()
    # Break if everything went fine
    break()
  endif()
endforeach()

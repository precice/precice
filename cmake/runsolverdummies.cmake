if(NOT EXISTS ${WRAPPER})
  message(FATAL_ERROR "CMake was unable to locate the execution wrapper script at ${WRAPPER}!")
endif()

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
  COMMAND ${WRAPPER} ${DUMMY_A} ${DUMMY_B} ${DUMMY_CONFIG}
  WORKING_DIRECTORY ${DUMMY_RUN_DIR}
  RESULT_VARIABLE DUMMY_RESULT
  )

# Check the return codes/statuses of the solvers
if(NOT (DUMMY_RESULT EQUAL 0))
  # Fail in case we encounter another error code/condition other than 0
  message(FATAL_ERROR "An error occured running the solverdummies! Return code : \"${DUMMY_RESULT}\"")
endif()

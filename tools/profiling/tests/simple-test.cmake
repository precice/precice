#
# Script to test event-precice merge trace and export given a folder
#
# Inputs:
#   FOLDER - the folder argument passed to merge
#   EVENTS_SCRIPT - the path to the precice-events script

find_package(Python3 COMPONENTS Interpreter REQUIRED)

if(NOT FOLDER OR NOT EVENTS_SCRIPT)
  message(FATAL_ERROR "Missing arguments")
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import pandas"
                RESULTS_VARIABLE PYTHON_NO_PANDAS)

message(STATUS "Removing previous files")
file(REMOVE events.json trace.json exports.csv)

message(STATUS "Testing: merge")
execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} merge ${FOLDER}
  COMMAND_ECHO STDOUT)
if(NOT EXISTS "events.json")
  message(FATAL_ERROR "No events.json file found")
else()
  if(NOT ${PYTHON_NO_PANDAS})
      message(STATUS "Testing: analyze")
      execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} analyze --output analyze.csv A events.json
      COMMAND_ECHO STDOUT)
      if(NOT EXISTS "analyze.csv")
        message(FATAL_ERROR "No analyze.csv file found")
      endif()
    endif()
endif()

message(STATUS "Testing: trace")
execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} trace
  COMMAND_ECHO STDOUT
  )
if(NOT EXISTS "trace.json")
  message(FATAL_ERROR "No trace.json file found")
endif()

message(STATUS "Testing: export")
execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} export
  COMMAND_ECHO STDOUT)
if(NOT EXISTS "events.csv")
  message(FATAL_ERROR "No events.csv file found")
endif()

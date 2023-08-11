#
# Script to test event-precice merge trace and export given a folder
#
# Reads the content of ${TEST_FOLDER}/.test to deduce which participants to analyze
#
# Inputs:
#   TEST_FOLDER - the folder of the tests
#   FOLDER_ARG - the folder argument passed to merge
#   EVENTS_SCRIPT - the path to the precice-events script

# Force using a virtual env if it is set in the environment
if(DEFINED ENV{VIRTUAL_ENV})
  set(Python3_FIND_VIRTUALENV ONLY)
endif()
find_package(Python3 COMPONENTS Interpreter REQUIRED)

if(NOT FOLDER_ARG OR NOT TEST_FOLDER OR NOT EVENTS_SCRIPT)
  message(FATAL_ERROR "Missing arguments")
endif()

message(STATUS "Removing previous files")
file(REMOVE events.json trace.json exports.csv)

message(STATUS "Testing: merge")
execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} merge ${FOLDER_ARG}
  COMMAND_ECHO STDOUT)
if(NOT EXISTS "events.json")
  message(FATAL_ERROR "No events.json file found")
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

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import pandas"
                RESULTS_VARIABLE PYTHON_NO_PANDAS)
if(NOT ${PYTHON_NO_PANDAS})
  file(STRINGS "${TEST_FOLDER}/.test" _participants)
  foreach(_participant IN LISTS _participants)
    message(STATUS "Testing: analyze ${_participant}")
    set(_outfile "analyze-${_participant}.csv")
    execute_process(COMMAND ${Python3_EXECUTABLE} ${EVENTS_SCRIPT} analyze --output ${_outfile} ${_participant} events.json
    COMMAND_ECHO STDOUT)
    if(NOT EXISTS "${_outfile}")
      message(FATAL_ERROR "No ${_outfile} file found")
    endif()
  endforeach()
endif()

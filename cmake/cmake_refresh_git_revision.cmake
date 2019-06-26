#
# This script uses git to extract the current status of the source control.
# It then configure the file version.cpp.in using this info.
#
# Required Arguments:
#   SRC - The full path to the configuration input (version.cpp.in)
#   DST - The full path to the configuration output (version.cpp)
#
# Configured Variables:
#   preCICE_COMMIT - the result of git describe --tags --broken
#

find_package(Git REQUIRED)
set(preCICE_COMMIT "no-info")
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --tags --broken
  RESULT_VARIABLE PRECICE_COMMIT_RET
  OUTPUT_VARIABLE preCICE_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
configure_file(${SRC} ${DST} @ONLY)

if("${PRECICE_COMMIT_RET}" EQUAL "0")
  message(STATUS "Revision status: ${preCICE_COMMIT}")
else()
  message(STATUS "Revision status: Detection failed")
endif()

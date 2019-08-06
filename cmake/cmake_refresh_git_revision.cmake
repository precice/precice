#
# This script uses git to extract the current status of the source control.
# It then configure the file version.cpp.in using this info.
#
# Required Arguments:
#   SRC - The full path to the configuration input (version.cpp.in)
#   DST - The full path to the configuration output (version.cpp)
#
# Configured Variables:
#   preCICE_REVISION - the result of git describe --tags --broken
#

find_package(Git REQUIRED)
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
  RESULT_VARIABLE PRECICE_REVISION_RET
  OUTPUT_VARIABLE PRECICE_REVISION_OUT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  )

set(preCICE_REVISION "no-info [Git failed/Not a repository]")
if("${PRECICE_REVISION_RET}" EQUAL "0")
  set(preCICE_REVISION "${PRECICE_REVISION_OUT}")
  message(STATUS "Revision status: ${preCICE_REVISION}")
else()
  message(STATUS "Revision status: Detection failed")
endif()
configure_file( ${SRC} ${DST} @ONLY)

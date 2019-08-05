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
set(preCICE_REVISION "no-info [Git failed/Not a repository]")
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --tags --broken
  RESULT_VARIABLE PRECICE_REVISION_RET
  OUTPUT_VARIABLE preCICE_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
configure_file(${SRC} ${DST} @ONLY)

if("${PRECICE_REVISION_RET}" EQUAL "0")
  message(STATUS "Revision status: ${preCICE_REVISION}")
else()
  message(STATUS "Revision status: Detection failed")
endif()

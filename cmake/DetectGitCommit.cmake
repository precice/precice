find_package(Git QUIET)
set(preCICE_COMMIT "no-info")
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --broken
    WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
    RESULT_VARIABLE PRECICE_COMMIT_RET
    OUTPUT_VARIABLE preCICE_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if("${PRECICE_COMMIT_RET}" EQUAL "0")
    message(STATUS "Revision status: ${preCICE_COMMIT}")
  else()
    message(WARNING "Revision status: Detection failed!")
  endif()
else(GIT_FOUND)
  message(STATUS "Revision status: No source version control detected.")
endif(GIT_FOUND)

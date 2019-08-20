#
# Git Revision Info
#
# This script integrates the git revision info into the preCICE library.
#
# If git is discoverable on the machine, then this script will register an
# additional target which is added as dependency to the precice target.
# The new target launches the script cmake_refresh_git_revision.cmake which
# automatically generates and updates the version.cpp file.
#
# If git was not found, this script simply configures a dummy versions.cpp file.
#

find_package(Git QUIET)
if(GIT_FOUND)
  add_custom_target(GitRevision
    COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR}/src/precice/impl/versions.cpp.in -DDST=${PROJECT_BINARY_DIR}/src/precice/impl/versions.cpp -DpreCICE_VERSION=${preCICE_VERSION} -DpreCICE_VERSION_INFORMATION=${preCICE_VERSION_INFORMATION} -P ${CMAKE_CURRENT_LIST_DIR}/cmake_refresh_git_revision.cmake
    WORKING_DIRECTORY "${preCICE_SOURCE_DIR}"
    BYPRODUCTS src/precice/impl/versions.cpp
    DEPENDS src/precice/impl/versions.cpp.in
    VERBATIM
    )
  add_dependencies(precice GitRevision)
else(GIT_FOUND)
  set(preCICE_REVISION "no-info [Git not found]")
  configure_file("${PROJECT_SOURCE_DIR}/src/precice/impl/versions.cpp.in" "${PROJECT_BINARY_DIR}/src/precice/impl/versions.cpp" @ONLY)
  message(STATUS "Revision status: No git installation detected.")
endif(GIT_FOUND)

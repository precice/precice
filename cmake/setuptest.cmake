#
# Setup script for the preCICE tests.
#
# The script wipes the given test directory.
#
# Arguments:
# DIR - the directory of the test
#
cmake_minimum_required(VERSION 3.10)

if(NOT DIR)
  message(FATAL_ERROR "You have to define a directory using the variable DIR !")
endif()

get_filename_component(_PATH ${DIR} ABSOLUTE)
file(REMOVE_RECURSE ${_PATH})
file(MAKE_DIRECTORY ${_PATH})

# This function writes an env file to be sourced by a shell.
# 
# Mandatory Arguments:
#   FILE      - The file to write
#   PREFIX    - The prefix for the options below
#
# Options:
#   LIBRARY   - The library dir, usually lib
#   RUNTIME   - The runtime dir, usually bin
#   INTERFACE - The interface dir, usually include
#
function(write_env_file)
  cmake_parse_arguments(PARSE_ARGV 0 WEF "" "FILE;LIBRARY;INTERFACE;RUNTIME" "")
  if(NOT WEF_FILE)
    message(FATAL_ERROR "The argument FILE is missing!")
  endif()
  if(NOT WEF_PREFIX)
    message(FATAL_ERROR "The argument PREFIX is missing!")
  endif()
  if(WEF_LIBRARY)
    file(APPEND ${WEF_FILE} "export LD_LIBRARY_PATH=\"${WEF_PREFIX}/${WEF_LIBRARY}:\$LD_LIBRARY_PATH\"\n")
    file(APPEND ${WEF_FILE} "export LIBRARY_PATH=\"${WEF_PREFIX}/${WEF_LIBRARY}:\$LIBRARY_PATH\"\n")
  endif()
  if(WEF_INTERFACE)
    file(APPEND ${WEF_FILE} "export CPATH=\"${WEF_PREFIX}/${WEF_INTERFACE}:\$CPATH\"\n")
  endif()
  if(WEF_RUNTIME)
    file(APPEND ${WEF_FILE} "export PATH=\"${WEF_PREFIX}/${WEF_RUNTIME}:\$PATH\"\n")
  endif()
endfunction(write_env_file)

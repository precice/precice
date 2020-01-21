# Prints a varaible VAR with a description DESC
function(print_variable VAR DESC)
  if(DEFINED ${VAR})
    if("${${VAR}}" STREQUAL "")
      message(STATUS "${DESC}: <empty>")
    else()
      message(STATUS "${DESC}: ${${VAR}}")
    endif()
  else()
    message(STATUS "${DESC}: <undefined>")
  endif()
endfunction(print_variable)


# Prints multiple variable-description tuples
# invoke: print_variables(VAR "VAR1;DESC1" "VAR2;DESC2")
function(print_variables)
  cmake_parse_arguments(PARSE_ARGV 0 PRINT_VARS "" "" "VARS")
  foreach(item ${PRINT_VARS_VARS})
    list(GET item 0 VAR)
    list(GET item 1 MSG)
    print_variable(${VAR} "${MSG}")
  endforeach(item)
endfunction(print_variables)


# Prints a single line
function(print_line)
  message(STATUS "--------------------")
endfunction(print_line)

function(print_empty)
  message("")
endfunction(print_empty)

# Prints a header in the log
function(print_section TITLE)
  message(STATUS "=== ${TITLE} ===")
endfunction(print_section)


# Prints a fixed set of general cmake configuration variables and additional variables
# invoke: to print the configuration only
#         print_configuration() 
# invoke: to print additional variables
#         print_configuration(ADDITIONAL "VAR1:DESC1" "VAR2:DESC2")
function(print_configuration)
  cmake_parse_arguments(PARSE_ARGV 0 PRINT_CONFIG "" "" "ADDITIONAL")
  print_section("CONFIGURATION")
  print_variables( VARS
    "PROJECT_VERSION;Library version to build"
    "CMAKE_BUILD_TYPE;Build configuration"
    "BUILD_SHARED_LIBS;Build shared libraries"
    "CMAKE_SYSTEM;Target system"
    "CMAKE_HOST_SYSTEM;Host system"
    "CMAKE_CXX_LIBRARY_ARCHITECTURE;Library architecture"
    "CMAKE_CXX_COMPILER;CXX compiler"
    "CMAKE_CXX_FLAGS;CXX compiler flags"
    "CMAKE_LINKER;CXX linker"
    "CMAKE_INSTALL_PREFIX;Install prefix"
    "CMAKE_SOURCE_DIR;Source directory"
    "CMAKE_BINARY_DIR;Binary directory"
    )
  if(PRINT_CONFIG_ADDITIONAL)
    print_variables(VARS ${PRINT_CONFIG_ADDITIONAL})
  endif()
endfunction(print_configuration)

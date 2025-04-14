# Prints a variable VAR with a description DESC
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
    "CMAKE_VERSION;CMake version"
    "PROJECT_VERSION;Library version to build"
  )

  if (CMAKE_CONFIGURATION_TYPES)
    print_variables( VARS "CMAKE_CONFIGURATION_TYPES;Build configurations")
  else()
    print_variables( VARS "CMAKE_BUILD_TYPE;Build configuration")
  endif()

  print_variables( VARS
    "BUILD_SHARED_LIBS;Build shared libraries"
    "CMAKE_SYSTEM;Target system"
    "CMAKE_HOST_SYSTEM;Host system"
    "CMAKE_INSTALL_PREFIX;Install prefix"
    "PROJECT_SOURCE_DIR;Source directory"
    "PROJECT_BINARY_DIR;Binary directory"
    "CMAKE_CXX_LIBRARY_ARCHITECTURE;Library architecture"
    "CMAKE_CXX_COMPILER;CXX compiler"
    "CMAKE_CXX_COMPILER_ID;CXX compiler ID"
    "CMAKE_CXX_COMPILER_VERSION;CXX compiler version"
    "CMAKE_CXX_FLAGS;CXX compiler flags"
    )

  if (CMAKE_CONFIGURATION_TYPES)
    foreach(type IN LISTS CMAKE_CONFIGURATION_TYPES)
      string(TOUPPER "${type}" _upper_build_type)
      print_variables( VARS "CMAKE_CXX_FLAGS_${_upper_build_type};CXX ${type} compiler flags")
    endforeach()
  else()
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _upper_build_type)
    print_variables( VARS "CMAKE_CXX_FLAGS_${_upper_build_type};CXX ${CMAKE_BUILD_TYPE} compiler flags")
  endif()

  if(CMAKE_VERSION VERSION_LESS 3.29)
    print_variables(VARS "CMAKE_LINKER;CXX linker")
  else()
    print_variables(VARS
    "CMAKE_CXX_COMPILER_LINKER;CXX linker"
    "CMAKE_CXX_COMPILER_LINKER_ID;CXX linker ID"
    "CMAKE_CXX_COMPILER_LINKER_VERSION;CXX linker version"
    )
  endif()
  print_variables(VARS
    "CMAKE_SHARED_LINKER_FLAGS;Shared linker flags"
    "CMAKE_EXE_LINKER_FLAGS;Executable linker flags"
  )
  if(PRINT_CONFIG_ADDITIONAL)
    print_variables(VARS ${PRINT_CONFIG_ADDITIONAL})
  endif()
endfunction(print_configuration)

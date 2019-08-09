#
# This file detects compilers and compiler flags set via the environment
# variables CPP and CPPFLAGS. This is necessary to be compliant with section 3
# of the xSDK Community Installation Policies
#
# !!! NOTE !!!
#   The compiler has to be set before any language is set in CMake.
#   This file has to be included BEFORE the first project() call!
#

# 3. Compilers and Flags
# a. i. Use env CPP as compiler
if(DEFINED ENV{CPP})
  set(CMAKE_CXX_COMPILER "$ENV{CPP}" CACHE STRING "The CXX compiler.")
endif()

# a. ii. Use env CPPFLAGS as compiler flags
if(DEFINED ENV{CPPFLAGS})
  if (CMAKE_CXX_FLAGS)
    list(APPEND CMAKE_CXX_FLAGS "$ENV{CPPFLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "$ENV{CPPFLAGS}" CACHE STRING "The cxx compiler flags.")
  endif()
endif()

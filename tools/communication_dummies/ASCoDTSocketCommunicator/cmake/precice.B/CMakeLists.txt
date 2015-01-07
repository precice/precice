# General {{{
# ==============================================================================
# Require minimum version of CMake.
cmake_minimum_required(VERSION 2.8.6)
# ------------------------------------------------------------------------------
# Disable in-source builds.
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()
# ------------------------------------------------------------------------------
# Declare project.
project("precice.B" C CXX)
# ------------------------------------------------------------------------------
# Issue deprecation errors for macros and functions.
set(CMAKE_ERROR_DEPRECATED ON)
# ------------------------------------------------------------------------------
# Be verbose and colorful.
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE   ON)
# ------------------------------------------------------------------------------
# Include the currently being processed both source and build directories.
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# ------------------------------------------------------------------------------
# Remove "lib" prefix for shared libraries (DLLs) on Windows.
if(WIN32)
  set(CMAKE_SHARED_LIBRARY_PREFIX "")
  set(CMAKE_SHARED_MODULE_PREFIX  "")
endif()
# ------------------------------------------------------------------------------
# Install to the current project by default.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../"
                           CACHE
                           PATH
                           ""
                           FORCE)
endif()
# ==============================================================================
# }}} General

# Packages {{{
# ==============================================================================
# JNI {{{
# ==============================================================================
find_package(JNI REQUIRED)

include_directories(${JNI_INCLUDE_DIRS})
# ==============================================================================
# }}} JNI

# Threads {{{
# ==============================================================================
set(CMAKE_THREAD_PREFER_PTHREAD ON)

find_package(Threads REQUIRED)

if(NOT CMAKE_USE_PTHREADS_INIT)
  message(FATAL_ERROR "POSIX threads not found.")
endif()

add_definitions(-DTHREADSAFE -D_REENTRANT)

set(Threads_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
# ==============================================================================
# }}} Threads

# Sockets {{{
# ==============================================================================
if(WIN32)
  set(Sockets_LIBRARIES ws2_32)
else()
  set(Sockets_LIBRARIES)
endif()
# ==============================================================================
# }}} Sockets

# MPI {{{
# ==============================================================================
if(PARALLEL)
  add_definitions(-DParallel)

  find_package(MPI REQUIRED)

  include_directories(${MPI_CXX_INCLUDE_PATH})

  set(MPI_LIBRARIES ${MPI_CXX_LIBRARIES})
endif()
# ==============================================================================
# }}} MPI
# ==============================================================================
# }}} Packages

# Compilers {{{
# ==============================================================================
if(CMAKE_COMPILER_IS_GNUCC)
  set(COMPILER_IS_GCC TRUE)
  set(LINKER_IS_LD    TRUE)
endif()
# ------------------------------------------------------------------------------
if(CMAKE_COMPILER_IS_GNUCXX)
  set(COMPILER_IS_G++ TRUE)
  set(LINKER_IS_LD    TRUE)
endif()
# ==============================================================================
# }}} Compilers

# Flags {{{
# ==============================================================================
if(COMPILER_IS_GCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fomit-frame-pointer")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -funroll-loops")
endif()
# ------------------------------------------------------------------------------
if(COMPILER_IS_G++)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fomit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")
endif()
# ------------------------------------------------------------------------------
if(LINKER_IS_LD)
  if(WIN32)
    if(COMPILER_IS_G++)
      set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -static-libstdc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static-libstdc++")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++")
    endif()

    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -static-libgcc")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static-libgcc")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc")

    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -static")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -static")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static")
  endif()
endif()
# ==============================================================================
# }}} Flags

# Targets {{{
# ==============================================================================
include_directories("../../src")
include_directories("../../components/c++")

include("sources.cmake")
include("targets.cmake")

target_link_libraries("precice.B"
                      ${JNI_LIBRARIES}
                      ${Threads_LIBRARIES}
                      ${Sockets_LIBRARIES}
                      ${MPI_LIBRARIES})

install(TARGETS "precice.B"
        ARCHIVE DESTINATION "lib"
        LIBRARY DESTINATION "lib"
        RUNTIME DESTINATION "bin")
# ==============================================================================
# }}} Targets

# Build the library if that hasn't happened yet

if(PRECICE_LIBBACKTRACE_PREPARED)
  message(STATUS "libbacktrace already build")
else()
  message(STATUS "Preparing libbacktrace")

  message(STATUS "Extracting sources")
  execute_process(
      COMMAND tar -xzf ${CMAKE_CURRENT_SOURCE_DIR}/libbacktrace.tar.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )

  set(_libbacktrace_sources ${CMAKE_CURRENT_BINARY_DIR}/libbacktrace-master)

  set(_step autoconf)
  message(STATUS "Running ${_step}")
  execute_process(
    COMMAND ${_step}
    WORKING_DIRECTORY ${_libbacktrace_sources}
    RESULTS_VARIABLE _r
    )
  message(STATUS "Step finished with ${_r}")

  message(STATUS "Running ./configure")
  execute_process(
    COMMAND ./configure --enable-shared=yes
    WORKING_DIRECTORY ${_libbacktrace_sources}
    RESULTS_VARIABLE _r
    )
  message(STATUS "Step finished with ${_r}")

  set(_step make)
  message(STATUS "Running ${_step}")
  execute_process(
    COMMAND ${_step}
    WORKING_DIRECTORY ${_libbacktrace_sources}
    RESULTS_VARIABLE _r
    )
  message(STATUS "Step finished with ${_r}")
  unset(_r)
  unset(step)

  message(STATUS "Copying files")
  file(COPY
    ${_libbacktrace_sources}/backtrace.h
    ${_libbacktrace_sources}/.libs/libbacktrace.so
    ${_libbacktrace_sources}/.libs/libbacktrace.so.0
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FOLLOW_SYMLINK_CHAIN)

  unset(_libbacktrace_sources)

  message(STATUS "libbacktrace built successfully")
  set(PRECICE_LIBBACKTRACE_PREPARED YES CACHE BOOL "Has the libbacktrace been build successfully?" FORCE)
  mark_as_advanced(PRECICE_LIBBACKTRACE_PREPARED)
endif()

# Create the library target
add_library(internal::libbacktrace SHARED IMPORTED GLOBAL)
set_target_properties(internal::libbacktrace PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/libbacktrace.so)
target_include_directories(internal::libbacktrace INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

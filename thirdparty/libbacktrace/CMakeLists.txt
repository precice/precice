# Build the library if that hasn't happened yet

if(PRECICE_LIBBACKTRACE_PREPARED)
  message(STATUS "libbacktrace already build")
else()
  message(STATUS "Preparing libbacktrace")

  message(STATUS "Extracting sources")
  file(ARCHIVE_EXTRACT
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/libbacktrace.zip
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )

  set(_libbacktrace_sources ${CMAKE_CURRENT_BINARY_DIR}/libbacktrace-master)

  foreach(_step autoconf ./configure make)
    message(STATUS "Running ${_step}")
    execute_process(
      COMMAND ${_step}
      WORKING_DIRECTORY ${_libbacktrace_sources}
      COMMAND_ERROR_IS_FATAL ANY
      )
  endforeach(_step)

  message(STATUS "Copying files")
  file(COPY
    ${_libbacktrace_sources}/backtrace.h
    ${_libbacktrace_sources}/.libs/libbacktrace.a
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  unset(_libbacktrace_sources)

  message(STATUS "libbacktrace built successfully")
  set(PRECICE_LIBBACKTRACE_PREPARED YES CACHE BOOL "Has the libbacktrace been build successfully?" FORCE)
  mark_as_advanced(PRECICE_LIBBACKTRACE_PREPARED)
endif()

# Create the library target
add_library(internal::libbacktrace STATIC IMPORTED GLOBAL)
set_target_properties(internal::libbacktrace PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/libbacktrace.a)
target_include_directories(internal::libbacktrace INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

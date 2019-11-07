#
# This function resolves a list of libraries based on a list of directories.
# Use this function to resolve the paths of a list of library names using a list of library directories.
# This is especially useful to handle the result of PkgConfig
#
function(resolve_library_locations PREFIX LIBS DIRS)  
  unset(_libs)
  foreach(lib ${LIBS})
    find_library(${PREFIX}_LIB_${lib}
      NAMES ${lib}
      PATHS ${DIRS}
      NO_DEFAULT_PATH
      )
    list(APPEND _libs ${${PREFIX}_LIB_${lib}})
  endforeach(lib)
  set(${PREFIX}_LIBRARIES_RESOLVED "${_libs}" PARENT_SCOPE)
endfunction(resolve_library_locations)


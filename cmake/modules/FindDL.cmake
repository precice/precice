# This module defines the following variables:
#
# DL_FOUND
#   True if libdl found.
#
# DL_LIBRARIES
#   The library libdl
#
# This module exports the following IMPORTED target:
#
# DL::DL
#

find_library(DL_LIBRARY
    NAMES dl libdl
    HINTS ${DL_ROOT} $ENV{DL_ROOT})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DL
    REQUIRED_VARS DL_LIBRARY)

if(DL_FOUND)
    set(DL_LIBRARIES ${DL_LIBRARY})

    if(NOT TARGET DL::DL)
        add_library(DL::DL UNKNOWN IMPORTED)
        set_target_properties(DL::DL PROPERTIES
            INTERFACE_LINK_LIBRARIES "${DL_LIBRARY}"
            )
        set_property(TARGET DL::DL APPEND PROPERTY IMPORTED_LOCATION "${DL_LIBRARY}")
    endif()
endif()

mark_as_advanced(DL_LIBRARIES)

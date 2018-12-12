# Cmake Finder module for PETSc.
#
# This will define a target PETSc which should include all necessary properties.
# Furthermore the variables PETSC_VERSION_MAJOR and PETSC_VERSION_MINOR will be set to their 
# respective values.
cmake_policy(VERSION 3.3)

find_package(PkgConfig REQUIRED)
set(ENV{PKG_CONFIG_PATH} "$ENV{CRAY_PETSC_PREFIX_DIR}/lib/pkgconfig:$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_search_module(PETSC craypetsc_real PETSc)
string(REPLACE "." ";" VERSION_LIST ${PETSC_VERSION})
list(GET VERSION_LIST 0 PETSC_VERSION_MAJOR)
list(GET VERSION_LIST 1 PETSC_VERSION_MINOR)
list(GET VERSION_LIST 2 PETSC_VERSION_SUBMINOR)
set (PETSC_VERSION_MAJOR ${PETSC_VERSION_MAJOR} CACHE STRING "Petsc Major version")
set (PETSC_VERSION_MINOR ${PETSC_VERSION_MINOR} CACHE STRING "Petsc Minor version")
add_library(PETSc INTERFACE IMPORTED)

list(LENGTH PETSC_LIBRARY_DIRS dir_count)
if(NOT ${dir_count} EQUAL 1)
    message(FATAL_ERROR "PETSc has more than one LIBRARY_DIR. This FindPETSc module is broken now.")
endif()

list(APPEND PETSC_LIBS_ABSOLUTE "")
foreach(libname ${PETSC_LIBRARIES})
    list(APPEND PETSC_LIBS_ABSOLUTE ${PETSC_LIBRARY_DIRS}/lib${libname}.so) 
endforeach()
set_property(TARGET PETSc PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
set_property(TARGET PETSc PROPERTY INTERFACE_LINK_LIBRARIES ${PETSC_LIBS_ABSOLUTE})

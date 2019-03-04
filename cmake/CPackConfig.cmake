#
# Packaging Settings for CPack
#
# Variable References:
#  General https://cmake.org/cmake/help/latest/module/CPack.html#module:CPack
#  Debian Packages https://cmake.org/cmake/help/latest/cpack_gen/deb.html#cpack_gen:CPack%20DEB%20Generator
#

# Install doc files
install(FILES tools/packaging/debian/copyright
  DESTINATION share/doc/libprecice${preCICE_VERSION}
  )

# Detect the system name
if(WIN32)
  set(CPACK_SYSTEM_NAME "win32")
elseif(WIN64)
  set(CPACK_SYSTEM_NAME "win64")
else()
  # Try to detect the codename of the distro using lsb_release
  find_program(LSB_RELEASE_EXE lsb_release)
  if(LSB_RELEASE_EXE)
    execute_process(COMMAND ${LSB_RELEASE_EXE} -cs
      OUTPUT_VARIABLE DISTRO_CODENAME
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    set(CPACK_SYSTEM_NAME "${DISTRO_CODENAME}")
  else()
    # Use the target system name of cmake as fallback
    set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
  endif()
endif()

# General
set(CPACK_PACKAGE_NAME "libprecice${preCICE_VERSION}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_SYSTEM_NAME}")
set(CPACK_PACKAGE_VENDOR "precice.org")
set(CPACK_PACKAGE_CONTACT "The precice developers <precice@mailman.informatik.uni-stuttgart.de>")
set(CPACK_PACKAGE_MAINTAINER "The precice developers <precice@mailman.informatik.uni-stuttgart.de>")
set(CPACK_PACKAGE_DESCRIPTION "preCICE (Precise Code Interaction Coupling Environment) is a coupling library for partitioned multi-physics simulations, including, but not restricted to fluid-structure interaction and conjugate heat transfer simulations. Partitioned means that preCICE couples existing programs (solvers) capable of simulating a subpart of the complete physics involved in a simulation. This allows for the high flexibility that is needed to keep a decent time-to-solution for complex multi-physics scenarios.")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Precise Code Interaction Coupling Environment")
set(CPACK_PACKAGE_EXECUTABLES "testprecice;binprecice")
set(CPACK_PACKAGE_HOMEPAGE_URL "www.precice.org")
#set(CPACK_PACKAGE_ICON "")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_RESOURCE_FILE_LICENSE "${preCICE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README  "${preCICE_SOURCE_DIR}/tools/packaging/README.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${preCICE_SOURCE_DIR}/tools/packaging/WELCOME.txt")
set(CPACK_MONOLITHIC_INSTALL TRUE)
set(CPACK_STRIP_FILES TRUE)
set(CPACK_GENERATOR "TGZ")
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr")
  list(APPEND CPACK_GENERATOR "DEB")
else()
  message("Debian package generator disabled: Install prefix is not \"/usr\"")
endif()

#set(CPACK_SOURCE_PACKAGE_FILE_NAME "")
set(CPACK_SOURCE_GENERATOR ${CPACK_GENERATOR})
set(CPACK_SOURCE_IGNORE_FILES
  "/build/"
  "/.git/"
  ".gitignore"
  )

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, petsc-dev (>= 3.6), libboost-dev (>= 1.65), libboost-log-dev (>= 1.65), libboost-thread-dev (>= 1.65), libboost-system-dev (>= 1.65), libboost-filesystem-dev (>= 1.65), libboost-program-options-dev (>= 1.65), libboost-test-dev (>= 1.65), libeigen3-dev, libxml2-dev, python-dev, python-numpy")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Precise Code Interaction Coupling Environment\n\
 preCICE (Precise Code Interaction Coupling Environment) is a coupling library\n\
 for partitioned multi-physics simulations, including, but not restricted to\n\
 fluid-structure interaction and conjugate heat transfer simulations.\n\
 Partitioned means that preCICE couples existing programs (solvers) capable of\n\
 simulating a subpart of the complete physics involved in a simulation.\n\
 This allows for the high flexibility that is needed to keep a decent\n\
 time-to-solution for complex multi-physics scenarios.\
")
set(CPACK_DEBIAN_PACKAGE_CONTROL_STRUCT_PERMISSION TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")

file(WRITE "${PRECICE_PACKAGING_DIR}/lintian-override" "${CPACK_PACKAGE_NAME} binary: non-dev-pkg-with-shlib-symlink")
install(FILES "${PRECICE_PACKAGING_DIR}/lintian-override" 
  DESTINATION share/lintian/overrides
  RENAME ${CPACK_PACKAGE_NAME}
  )

# Configure a pkg-config file for the debian package
configure_file(
  "${PROJECT_SOURCE_DIR}/tools/packaging/debian/precice.pc.in"
  "${PRECICE_PACKAGING_DIR}/pkgconfig/libprecice.pc"
  @ONLY
  )
install(DIRECTORY "${PRECICE_PACKAGING_DIR}/pkgconfig" 
  DESTINATION lib
)


include(CPack)

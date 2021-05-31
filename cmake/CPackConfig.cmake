#
# Packaging Settings for CPack
#
# Variable References:
#  General https://cmake.org/cmake/help/latest/module/CPack.html#module:CPack
#  Debian Packages https://cmake.org/cmake/help/latest/cpack_gen/deb.html#cpack_gen:CPack%20DEB%20Generator
#

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
set(CPACK_PACKAGE_NAME "libprecice${preCICE_SOVERSION}")
set(CPACK_PACKAGE_VERSION "${preCICE_VERSION}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CPACK_SYSTEM_NAME}")
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
set(CPACK_RESOURCE_FILE_README  "${preCICE_SOURCE_DIR}/tools/releasing/packaging/README.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${preCICE_SOURCE_DIR}/tools/releasing/packaging/WELCOME.txt")
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

# Build dependecy set
unset(CPACK_DEBIAN_PACKAGE_DEPENDS)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libboost-dev (>= 1.65), libboost-log-dev (>= 1.65), libboost-thread-dev (>= 1.65), libboost-system-dev (>= 1.65), libboost-filesystem-dev (>= 1.65), libboost-program-options-dev (>= 1.65), libboost-test-dev (>= 1.65), libxml2")
if(PRECICE_PythonActions)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python3-dev, python3-numpy")
endif()
if(PRECICE_MPICommunication)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, mpi-default-dev")
endif()
if(PRECICE_PETScMapping)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev (>= 3.6)")
endif()

set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "\
 preCICE is a coupling library for partitioned multi-physics simulations,\n\
 including, but not restricted to fluid-structure interaction and\n\
 conjugate heat transfer simulations.\n\
 Partitioned means that preCICE couples existing programs (solvers) capable of\n\
 simulating a subpart of the complete physics involved in a simulation.\n\
 This allows for the high flexibility that is needed to keep a decent\n\
 time-to-solution for complex multi-physics scenarios.\
")
set(CPACK_DEBIAN_PACKAGE_CONTROL_STRUCT_PERMISSION TRUE)
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${preCICE_SOURCE_DIR}/tools/releasing/packaging/debian/triggers")
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")

# Install doc files
install(FILES tools/releasing/packaging/debian/copyright
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/${CPACK_PACKAGE_NAME}
  )

# Install lintian override
file(WRITE "${PRECICE_PACKAGING_DIR}/lintian-override" "${CPACK_PACKAGE_NAME} binary: non-dev-pkg-with-shlib-symlink")
install(FILES "${PRECICE_PACKAGING_DIR}/lintian-override" 
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/lintian/overrides
  RENAME ${CPACK_PACKAGE_NAME}
  )

# Compress and install the debian changelog
find_program(GZIP_EXE gzip DOC "The gzip executable")
if(GZIP_EXE)
  # Process the changelog for debian package
  message(STATUS "Compressing changelog")
  file(COPY tools/releasing/packaging/debian/changelog DESTINATION ${PRECICE_PACKAGING_DIR})
  execute_process(COMMAND "${GZIP_EXE}" "-9nf" "${PRECICE_PACKAGING_DIR}/changelog")

  # Install compressed changelog
  install(FILES ${PRECICE_PACKAGING_DIR}/changelog.gz
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/${CPACK_PACKAGE_NAME}
    )
else()
  message(WARNING "Installing uncompressed changelog")
  # Install uncompressed changelog
  install(FILES tools/releasing/packaging/debian/changelog
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/${CPACK_PACKAGE_NAME}
    )
endif()

file(REMOVE CPackConfig.cmake CPackSourceConfig.cmake)
include(CPack)

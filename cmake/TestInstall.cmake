#
# Call the test_install target to test the installation of the library
#
# This configures, builds and tests the c++ solverdummy in the TestInstall directory
#
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TestInstall)
add_custom_target(test_install
  COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/TestInstall/CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/TestInstall/CMakeFiles
  COMMAND ${CMAKE_COMMAND} -Dprecice_DIR=${CMAKE_INSTALL_PREFIX}/lib/cmake/precice ${preCICE_SOURCE_DIR}/examples/solverdummies/cpp
  COMMAND ${CMAKE_COMMAND} --build . --target all
  COMMAND ${CMAKE_CTEST_COMMAND} -V
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TestInstall
  COMMENT "Testing the installation using the C++ solverdummy"
  VERBATIM
  )

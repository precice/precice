# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${preCICE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${preCICE_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${preCICE_BINARY_DIR}/cmake_uninstall.cmake)
endif()

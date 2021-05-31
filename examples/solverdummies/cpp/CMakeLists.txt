cmake_minimum_required(VERSION 3.10.2)
project(SolverdummyCPP LANGUAGES CXX DESCRIPTION "preCICE CPP solver dummy")

find_package(precice REQUIRED CONFIG)

add_executable(solverdummy solverdummy.cpp)
set_target_properties(solverdummy PROPERTIES CXX_STANDARD 11)
target_link_libraries(solverdummy PRIVATE precice::precice)


enable_testing()
add_test(NAME precice.solverdummy
  COMMAND ${CMAKE_COMMAND}
  -DDUMMY_EXE=$<TARGET_FILE:solverdummy>
  -DDUMMY_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/../precice-config.xml
  -V
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )

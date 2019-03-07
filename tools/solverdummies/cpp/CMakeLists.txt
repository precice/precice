cmake_minimum_required(VERSION 3.9)
project(SolverdummyCPP LANGUAGES CXX DESCRIPTION "preCICE CPP solver dummy")

find_package(precice REQUIRED CONFIG)

add_executable(solverdummy solverdummy.cpp)
set_target_properties(solverdummy PROPERTIES CXX_STANDARD 11)
target_link_libraries(solverdummy PRIVATE precice::precice)

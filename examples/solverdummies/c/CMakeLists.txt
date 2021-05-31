cmake_minimum_required(VERSION 3.10.2)
project(SolverdummyC LANGUAGES C DESCRIPTION "preCICE C solver dummy")

find_package(precice REQUIRED CONFIG)

add_executable(solverdummy solverdummy.c)
target_link_libraries(solverdummy PRIVATE precice::precice)

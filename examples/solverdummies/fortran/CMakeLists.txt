cmake_minimum_required(VERSION 3.10.2)
project(SolverdummyFortran LANGUAGES C Fortran DESCRIPTION "preCICE Fortran solver dummy")

find_package(precice REQUIRED CONFIG)

add_executable(solverdummy solverdummy.f90)
target_link_libraries(solverdummy PRIVATE precice::precice)

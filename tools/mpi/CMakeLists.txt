cmake_minimum_required (VERSION 3.10.2)

project(MPIPortsTester VERSION 1.0.0 LANGUAGES CXX)

find_package(MPI REQUIRED)

add_executable(acceptor acceptor/main.cpp)
target_link_libraries(acceptor PRIVATE MPI::MPI_CXX)

add_executable(requester requester/main.cpp)
target_link_libraries(requester PRIVATE MPI::MPI_CXX)

set_property(TARGET precice PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/SolverInterfaceFortran.hpp
  APPEND)

target_sources(precice PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/src/SolverInterfaceFortran.cpp
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/SolverInterfaceFortran.hpp
  )

target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )

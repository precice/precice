set_property(TARGET precice PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFortran.hpp
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFASTEST.hpp
  APPEND)

target_sources(precice PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFortran.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFASTEST.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFortran.hpp
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceFASTEST.hpp
  )

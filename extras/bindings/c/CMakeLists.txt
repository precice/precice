set_property(TARGET precice PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceC.h
  APPEND)

target_sources(precice PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceC.h
  ${CMAKE_CURRENT_LIST_DIR}/SolverInterfaceC.cpp
  )

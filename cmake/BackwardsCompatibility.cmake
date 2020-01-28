# REMOVE THIS FILE AS SOON AS TESTS ARE FIXED

if(DEFINED MPI)
  message(WARNING "The MPI flag was deprecated, please use PRECICE_MPICommunication instead.")
  set(PRECICE_MPICommunication ${MPI} CACHE BOOL "" FORCE)
endif()

if(DEFINED PETSC)
  message(WARNING "The PETSC flag was deprecated, please use PRECICE_PETScMapping instead.")
  set(PRECICE_PETScMapping ${PETSC} CACHE BOOL "" FORCE)
endif()

if(DEFINED PYTHON)
  message(WARNING "The PYTHON flag was deprecated, please use PRECICE_PythonActions instead.")
  set(PRECICE_PythonActions ${PYTHON} CACHE BOOL "" FORCE)
endif()

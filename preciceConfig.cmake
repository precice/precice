include(CMakeFindDependencyMacro)

# Required
find_dependency(Boost 1.60.0)
find_dependency(Eigen3 3.2)
find_dependency(LibXml2)

# Configurable
find_dependency(MPI)
find_dependency(PETSc 3.6)
find_dependency(PythonLibs 2.7)
find_dependency(NumPy)

include("${CMAKE_CURRENT_LIST_DIR}/preciceTargets.cmake")

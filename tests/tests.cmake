#
# This file lists all integration test sources and test suites
#
target_sources(testprecice
    PRIVATE
    tests/parallel/GlobalRBFPartitioning.cpp
    tests/parallel/LocalRBFPartitioning.cpp
    tests/parallel/TestFinalize.cpp
    tests/parallel/lifecycle/ConstructAndExplicitFinalize.cpp
    tests/parallel/lifecycle/ConstructOnly.cpp
    tests/parallel/lifecycle/Full.cpp
    tests/parallel/lifecycle/ImplicitFinalize.cpp
    tests/serial/MultiCouplingFourSolvers1.cpp
    tests/serial/MultiCouplingFourSolvers2.cpp
    tests/serial/helpers.cpp
    tests/serial/helpers.hpp
    tests/serial/initialize-data/Explicit.cpp
    tests/serial/initialize-data/Implicit.cpp
    tests/serial/initialize-data/ReadMapping.cpp
    tests/serial/initialize-data/WriteMapping.cpp
    tests/serial/initialize-data/helpers.cpp
    tests/serial/initialize-data/helpers.hpp
    tests/serial/mesh-requirements/NearestNeighborA.cpp
    tests/serial/mesh-requirements/NearestNeighborB.cpp
    tests/serial/mesh-requirements/NearestProjection2DA.cpp
    tests/serial/mesh-requirements/NearestProjection2DB.cpp
    tests/serial/time/explicit/DoNothingWithSubcycling.cpp
    tests/serial/time/explicit/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/whitebox/TestConfigurationComsol.cpp
    tests/serial/whitebox/TestConfigurationPeano.cpp
    tests/serial/whitebox/TestExplicitWithDataScaling.cpp
    )

# Contains the list of integration test suites
set(PRECICE_TEST_SUITES Parallel Serial)

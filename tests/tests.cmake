#
# This file lists all integration test sources and test suites
#
target_sources(testprecice
    PRIVATE
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
    tests/serial/time/explicit/DoNothingWithSubcycling.cpp
    tests/serial/time/explicit/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/ReadWriteScalarDataWithSubcycling.cpp
    )

# Contains the list of integration test suites
set(PRECICE_TEST_SUITES Serial)

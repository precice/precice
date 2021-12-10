#
# This file lists all tests sources that will be compiled into the test executable
#
target_sources(testprecice
    PRIVATE
    tests/serial/MultiCouplingFourSolvers1.cpp
    tests/serial/MultiCouplingFourSolvers2.cpp
    tests/serial/helpers.cpp
    tests/serial/helpers.hpp
    tests/serial/initialize-data/helpers.cpp
    tests/serial/initialize-data/helpers.hpp
    tests/serial/initialize-data/testDataInitializationReadMapping.cpp
    tests/serial/initialize-data/testDataInitializationWriteMapping.cpp
    tests/serial/initialize-data/testExplicitWithDataInitialization.cpp
    tests/serial/initialize-data/testImplicitWithDataInitialization.cpp
    )

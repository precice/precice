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
    tests/serial/time/explicit/compositional/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/explicit/serial-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSamplingFirst.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSamplingFirstNoInit.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSamplingZero.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingFirst.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingZero.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithSubcycling.cpp
    )

# Contains the list of integration test suites
set(PRECICE_TEST_SUITES Serial)

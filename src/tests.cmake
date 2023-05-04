#
# This file lists all tests sources that will be compiled into the test executable
#
target_sources(testprecice
    PRIVATE
    src/acceleration/test/AccelerationIntraCommTest.cpp
    src/acceleration/test/AccelerationSerialTest.cpp
    src/acceleration/test/ParallelMatrixOperationsTest.cpp
    src/acceleration/test/PreconditionerTest.cpp
    src/acceleration/test/QRFactorizationTest.cpp
    src/action/tests/PythonActionTest.cpp
    src/action/tests/ScaleActionTest.cpp
    src/action/tests/SummationActionTest.cpp
    src/com/tests/CommunicateBoundingBoxTest.cpp
    src/com/tests/CommunicateMeshTest.cpp
    src/com/tests/GenericTestFunctions.hpp
    src/com/tests/MPIDirectCommunicationTest.cpp
    src/com/tests/MPIPortsCommunicationTest.cpp
    src/com/tests/MPISinglePortsCommunicationTest.cpp
    src/com/tests/SocketCommunicationTest.cpp
    src/cplscheme/tests/AbsoluteConvergenceMeasureTest.cpp
    src/cplscheme/tests/CompositionalCouplingSchemeTest.cpp
    src/cplscheme/tests/DummyCouplingScheme.cpp
    src/cplscheme/tests/DummyCouplingScheme.hpp
    src/cplscheme/tests/ExplicitCouplingSchemeTest.cpp
    src/cplscheme/tests/MinIterationConvergenceMeasureTest.cpp
    src/cplscheme/tests/ParallelImplicitCouplingSchemeTest.cpp
    src/cplscheme/tests/RelativeConvergenceMeasureTest.cpp
    src/cplscheme/tests/ResidualRelativeConvergenceMeasureTest.cpp
    src/cplscheme/tests/SerialImplicitCouplingSchemeTest.cpp
    src/io/tests/ExportCSVTest.cpp
    src/io/tests/ExportConfigurationTest.cpp
    src/io/tests/ExportVTKTest.cpp
    src/io/tests/ExportVTPTest.cpp
    src/io/tests/ExportVTUTest.cpp
    src/io/tests/TXTTableWriterTest.cpp
    src/io/tests/TXTWriterReaderTest.cpp
    src/m2n/tests/GatherScatterCommunicationTest.cpp
    src/m2n/tests/PointToPointCommunicationTest.cpp
    src/mapping/tests/LinearCellInterpolationMappingTest.cpp
    src/mapping/tests/MappingConfigurationTest.cpp
    src/mapping/tests/NearestNeighborGradientMappingTest.cpp
    src/mapping/tests/NearestNeighborMappingTest.cpp
    src/mapping/tests/NearestProjectionMappingTest.cpp
    src/mapping/tests/PartitionOfUnityClusteringTest.cpp
    src/mapping/tests/PartitionOfUnityMappingTest.cpp
    src/mapping/tests/PetRadialBasisFctMappingTest.cpp
    src/mapping/tests/PolationTest.cpp
    src/mapping/tests/RadialBasisFctMappingTest.cpp
    src/math/tests/BarycenterTest.cpp
    src/math/tests/DifferencesTest.cpp
    src/math/tests/GeometryTest.cpp
    src/math/tests/MathTest.cpp
    src/mesh/tests/BoundingBoxTest.cpp
    src/mesh/tests/DataConfigurationTest.cpp
    src/mesh/tests/EdgeTest.cpp
    src/mesh/tests/FilterTest.cpp
    src/mesh/tests/MeshTest.cpp
    src/mesh/tests/TetrahedronTest.cpp
    src/mesh/tests/TriangleTest.cpp
    src/mesh/tests/VertexTest.cpp
    src/partition/tests/ProvidedPartitionTest.cpp
    src/partition/tests/ReceivedPartitionTest.cpp
    src/partition/tests/fixtures.hpp
    src/precice/tests/DataContextTest.cpp
    src/precice/tests/ParallelTests.cpp
    src/precice/tests/SpanTests.cpp
    src/precice/tests/ToolingTests.cpp
    src/precice/tests/VersioningTests.cpp
    src/precice/tests/WatchIntegralTest.cpp
    src/precice/tests/WatchPointTest.cpp
    src/query/tests/RTreeAdapterTests.cpp
    src/query/tests/RTreeTests.cpp
    src/testing/DataContextFixture.cpp
    src/testing/DataContextFixture.hpp
    src/testing/GlobalFixtures.cpp
    src/testing/ParallelCouplingSchemeFixture.cpp
    src/testing/ParallelCouplingSchemeFixture.hpp
    src/testing/SerialCouplingSchemeFixture.cpp
    src/testing/SerialCouplingSchemeFixture.hpp
    src/testing/TestContext.cpp
    src/testing/TestContext.hpp
    src/testing/Testing.cpp
    src/testing/Testing.hpp
    src/testing/WaveformFixture.cpp
    src/testing/WaveformFixture.hpp
    src/testing/main.cpp
    src/testing/tests/ExampleTests.cpp
    src/time/tests/StorageTest.cpp
    src/time/tests/WaveformTest.cpp
    src/utils/tests/AlgorithmTest.cpp
    src/utils/tests/DimensionsTest.cpp
    src/utils/tests/EigenHelperFunctionsTest.cpp
    src/utils/tests/IntraCommTest.cpp
    src/utils/tests/ManageUniqueIDsTest.cpp
    src/utils/tests/MultiLockTest.cpp
    src/utils/tests/ParallelTest.cpp
    src/utils/tests/StatisticsTest.cpp
    src/utils/tests/StringTest.cpp
    src/xml/tests/ParserTest.cpp
    src/xml/tests/PrinterTest.cpp
    src/xml/tests/XMLTest.cpp
    )

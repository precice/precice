#
# This file lists all tests sources that will be compiled into the test executable
#
target_sources(testprecice
    PRIVATE
    src/math/tests/GeometryTest.cpp
    src/math/tests/DifferencesTest.cpp
    src/math/tests/BarycenterTest.cpp
    src/partition/tests/ProvidedPartitionTest.cpp
    src/partition/tests/ReceivedPartitionTest.cpp
    src/query/tests/FindVoxelContentTest.cpp
    src/query/tests/FindClosestVertexVisitorTest.cpp
    src/query/tests/FindClosestTest.cpp
    src/m2n/tests/GatherScatterCommunicationTest.cpp
    src/m2n/tests/PointToPointCommunicationTest.cpp
    src/cplscheme/tests/SerialImplicitCouplingSchemeTest.cpp
    src/cplscheme/tests/DummyCouplingScheme.hpp
    src/cplscheme/tests/ParallelImplicitCouplingSchemeTest.cpp
    src/cplscheme/tests/CompositionalCouplingSchemeTest.cpp
    src/cplscheme/tests/PreconditionerTest.cpp
    src/cplscheme/tests/PostProcessingMasterSlaveTest.cpp
    src/cplscheme/tests/QRFactorizationTest.cpp
    src/cplscheme/tests/MinIterationConvergenceMeasureTest.cpp
    src/cplscheme/tests/AbsoluteConvergenceMeasureTest.cpp
    src/cplscheme/tests/HierarchicalAitkenPostProcessingTest.cpp
    src/cplscheme/tests/DummyCouplingScheme.cpp
    src/cplscheme/tests/ExplicitCouplingSchemeTest.cpp
    src/cplscheme/tests/RelativeConvergenceMeasureTest.cpp
    src/cplscheme/tests/ParallelMatrixOperationsTest.cpp
    src/action/tests/ScaleActionTest.cpp
    src/action/tests/PythonActionTest.cpp
    src/action/tests/ModifyCoordinatesActionTest.cpp
    src/utils/tests/ParallelTest.cpp
    src/utils/tests/StringTest.cpp
    src/utils/tests/DimensionsTest.cpp
    src/utils/tests/ManageUniqueIDsTest.cpp
    src/utils/tests/PointerVectorTest.cpp
    src/io/tests/ExportConfigurationTest.cpp
    src/io/tests/ExportVTKXMLTest.cpp
    src/io/tests/TXTTableWriterTest.cpp
    src/io/tests/TXTWriterReaderTest.cpp
    src/io/tests/ExportVTKTest.cpp
    src/com/tests/CommunicateMeshTest.cpp
    src/com/tests/MPIDirectCommunicationTest.cpp
    src/com/tests/SocketCommunicationTest.cpp
    src/com/tests/CommunicateBoundingBoxTest.cpp
    src/com/tests/MPIPortsCommunicationTest.cpp
    src/com/tests/GenericTestFunctions.hpp
    src/testing/main.cpp
    src/testing/Testing.cpp
    src/testing/Fixtures.hpp
    src/testing/Testing.hpp
    src/testing/tests/ExampleTests.cpp
    src/mesh/tests/EdgeTest.cpp
    src/mesh/tests/RTreeTests.cpp
    src/mesh/tests/MergeTest.cpp
    src/mesh/tests/MeshTest.cpp
    src/mesh/tests/PropertyContainerTest.cpp
    src/mesh/tests/TriangleTest.cpp
    src/mesh/tests/VertexTest.cpp
    src/mesh/tests/QuadTest.cpp
    src/mesh/tests/GroupTest.cpp
    src/mesh/tests/DataConfigurationTest.cpp
    src/xml/tests/XMLTest.cpp
    src/xml/tests/ParserTest.cpp
    src/precice/tests/MeshHandleTest.cpp
    src/precice/tests/WatchPointTest.cpp
    src/precice/tests/ServerTests.cpp
    src/precice/tests/ParallelTests.cpp
    src/precice/tests/SerialTests.cpp
    src/mapping/tests/NearestNeighborMappingTest.cpp
    src/mapping/tests/MappingConfigurationTest.cpp
    src/mapping/tests/RadialBasisFctMappingTest.cpp
    src/mapping/tests/PetRadialBasisFctMappingTest.cpp
    src/mapping/tests/NearestProjectionMappingTest.cpp
    )

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <precice/impl/SolverInterfaceImpl.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(TestOneTriangle)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using precice::testing::equals;

  // Implement your test here.
  
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<precice::VertexID> vertexIDs;
  

  if (context.isNamed("SolverOne")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);

    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());
    
    auto edgeAB = interface.setMeshEdge(meshID, vertexIDs[0], vertexIDs[1]);
    auto edgeBC = interface.setMeshEdge(meshID, vertexIDs[1], vertexIDs[2]);
    auto edgeCA = interface.setMeshEdge(meshID, vertexIDs[2], vertexIDs[0]);

    interface.setMeshTriangle(meshID, edgeAB, edgeBC, edgeCA);

    BOOST_CHECK(interface.getMeshVertexSize(meshID) == 3);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 3);
    BOOST_REQUIRE(mesh.edges().size() == 3);
    BOOST_REQUIRE(mesh.triangles().size() == 1);
    
    // Initialize, write data, advance and finalize
    double dt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values{1.0, 3.0, 7.0};
    interface.writeBlockScalarData(dataID, 3, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();

  } else {
    auto meshID = interface.getMeshID("MeshTwo");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords{1. / 3, 1. / 3};
    vertexIDs.resize(coords.size() / 2);

    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    // Initialize, read data, advance and finalize. Check expected mapping
    double dt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    //Check expected VS read
    Eigen::VectorXd expected(1);
    Eigen::VectorXd readData(1);
    expected << 11.0 / 3;

    interface.readBlockScalarData(dataID, expected.size(), vertexIDs.data(), readData.data());

    BOOST_CHECK(equals(expected, readData));

    interface.finalize();

  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI

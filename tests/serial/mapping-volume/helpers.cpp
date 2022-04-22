#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "mesh/Utils.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"

void testMappingVolumeOneTriangle(const std::string configFile, const TestContext &context)
{
  using precice::testing::equals;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<precice::VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords{0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);

    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    BOOST_TEST(vertexIDs[0] != -1, "Vertex A is invalid");
    BOOST_TEST(vertexIDs[1] != -1, "Vertex B is invalid");
    BOOST_TEST(vertexIDs[2] != -1, "Vertex C is invalid");

    int edgeAB = interface.setMeshEdge(meshID, vertexIDs[0], vertexIDs[1]);
    int edgeBC = interface.setMeshEdge(meshID, vertexIDs[1], vertexIDs[2]);
    int edgeCA = interface.setMeshEdge(meshID, vertexIDs[2], vertexIDs[0]);

    BOOST_TEST(edgeAB != -1, "Edge AB is invalid");
    BOOST_TEST(edgeBC != -1, "Edge BC is invalid");
    BOOST_TEST(edgeCA != -1, "Edge CA is invalid");

    interface.setMeshTriangle(meshID, edgeAB, edgeBC, edgeCA);

    BOOST_CHECK(interface.getMeshVertexSize(meshID) == 3);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 3);
    BOOST_REQUIRE(mesh.edges().size() == 3);
    BOOST_REQUIRE(mesh.triangles().size() == 1);

    BOOST_TEST(equals(mesh.triangles()[0].getArea(), 0.5), "Triangle area must be 0.5");

    // Initialize, write data, advance and finalize
    double dt = interface.initialize();
    BOOST_TEST(!mesh.triangles().empty(), "Triangle must still be stored");
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values{1.0, 100.0, 10.0};
    interface.writeBlockScalarData(dataID, 3, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();

  } else {
    auto meshID = interface.getMeshID("MeshTwo");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords{1. / 3., 1. / 3.};
    vertexIDs.resize(coords.size() / 2);

    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    // Initialize, read data, advance and finalize. Check expected mapping
    double dt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    // If "read" mapping, check received mesh
    if (precice::testing::WhiteboxAccessor::impl(interface).hasMesh("MeshOne")) {
      auto &mesh = precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
      BOOST_REQUIRE(mesh.vertices().size() == 3);
      BOOST_REQUIRE(mesh.edges().size() == 3);
      BOOST_REQUIRE(mesh.triangles().size() == 1);
    }

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    //Check expected VS read
    Eigen::VectorXd expected(1);
    Eigen::VectorXd readData(1);
    expected << 111.0 / 3;

    interface.readBlockScalarData(dataID, expected.size(), vertexIDs.data(), readData.data());
    // Current problem: "No triangle found", despite correct triangle set!
    BOOST_CHECK(equals(expected, readData));

    interface.finalize();
  }
}

#endif
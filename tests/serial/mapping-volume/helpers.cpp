#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/precice.hpp"

void testMappingVolumeOneTriangle(const std::string configFile, const TestContext &context, bool read)
{
  using precice::testing::equals;

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<precice::VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords{0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    BOOST_TEST(vertexIDs[0] != -1, "Vertex A is invalid");
    BOOST_TEST(vertexIDs[1] != -1, "Vertex B is invalid");
    BOOST_TEST(vertexIDs[2] != -1, "Vertex C is invalid");

    participant.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2]);

    BOOST_CHECK(participant.getMeshVertexSize(meshName) == 3);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    BOOST_REQUIRE(mesh.nVertices() == 3);
    BOOST_REQUIRE(mesh.triangles().size() == 1);

    BOOST_TEST(equals(mesh.triangles()[0].getArea(), 0.5), "Triangle area must be 0.5");

    // Initialize, write data, advance and finalize
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(!mesh.triangles().empty(), "Triangle must still be stored");
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values{1.0, 100.0, 10.0};
    participant.writeData(meshName, dataName, vertexIDs, values);

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant must advance only once.");
    participant.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords{1. / 3., 1. / 3.};
    vertexIDs.resize(coords.size() / 2);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    // Initialize, read data, advance and finalize. Check expected mapping
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant must advance once.");

    // If "read" mapping, check received mesh
    if (read) {
      auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
      BOOST_REQUIRE(mesh.nVertices() == 3);
      BOOST_REQUIRE(mesh.edges().size() == 3);
      BOOST_REQUIRE(mesh.triangles().size() == 1);
    }

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(1);
    Eigen::VectorXd readData(1);
    // Expected value in the middle of the triangle is the average of inputs (1, 10, 100)
    expected << 111.0 / 3;
    dt = participant.getMaxTimeStepSize();
    participant.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));

    participant.finalize();
  }
}

void testMappingVolumeOneTriangleConservative(const std::string configFile, const TestContext &context)
{
  using precice::testing::equals;

  precice::Participant participant(context.name, context.config(), context.rank, context.size);
  // SolverOne defines a vertex and a conserved quantity (e.g. a force) on it.
  // SolverTwo defines a triangle and read the mapped quantity. We check it is spread correctly.

  std::vector<precice::VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords{0.3, 0.2};
    vertexIDs.resize(coords.size() / 2);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    BOOST_TEST(vertexIDs[0] != -1, "Vertex A is invalid");
    BOOST_CHECK(participant.getMeshVertexSize(meshName) == 1);

    // Initialize, write data, advance and finalize
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values{1.0};
    participant.writeData(meshName, dataName, vertexIDs, values);

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant must advance only once.");
    participant.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords{0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    participant.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2]);

    // Initialize, read data, advance and finalize. Check expected mapping
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant must advance once.");

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(3);
    Eigen::VectorXd readData(3);
    // For conservative load, each point takes a fraction of the load.
    // This fraction is the barycentric coordinate, and there the load is 1.
    // Input point is (0.3, 0.2) and barycentric coordinates are thus (0.5, 0.3, 0.2)
    expected << 0.5, 0.3, 0.2;

    dt = participant.getMaxTimeStepSize();
    participant.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));

    participant.finalize();
  }
}

void testMappingVolumeOneTetra(const std::string configFile, const TestContext &context, bool read)
{
  using precice::testing::equals;

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<precice::VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords{0.0, 0.0, 0.0,
                               1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 3);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    BOOST_TEST(vertexIDs[0] != -1, "Vertex A is invalid");
    BOOST_TEST(vertexIDs[1] != -1, "Vertex B is invalid");
    BOOST_TEST(vertexIDs[2] != -1, "Vertex C is invalid");
    BOOST_TEST(vertexIDs[3] != -1, "Vertex D is invalid");

    participant.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);

    BOOST_CHECK(participant.getMeshVertexSize(meshName) == 4);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    // setMeshTetrahedron currently adds underlying connectivity
    BOOST_REQUIRE(mesh.nVertices() == 4);
    BOOST_REQUIRE(mesh.edges().empty());
    BOOST_REQUIRE(mesh.triangles().empty());
    BOOST_REQUIRE(mesh.tetrahedra().size() == 1);

    BOOST_TEST(equals(mesh.tetrahedra()[0].getVolume(), 1.0 / 6), "Tetrahedron volume must be 1/6");

    // Initialize, write data, advance and finalize
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_REQUIRE(mesh.edges().size() == 6);
    BOOST_REQUIRE(mesh.triangles().size() == 4);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 1);
    BOOST_TEST(!mesh.tetrahedra().empty(), "Tetrahedra must still be stored");
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant must advance once.");

    // Send 1 + 5x - 3y + 9z
    std::vector<double> values{1.0, 6.0, -2.0, 8.0};
    participant.writeData(meshName, dataName, vertexIDs, values);

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant must advance only once.");
    participant.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords{0.25, 0.25, 0.25};
    vertexIDs.resize(coords.size() / 2);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    // Initialize, read data, advance and finalize. Check expected mapping
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant must advance once.");

    // If "read" mapping, check received mesh, including connectivity
    if (read) {
      auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
      BOOST_CHECK(mesh.nVertices() == 4);
      BOOST_CHECK(mesh.edges().size() == 6);
      BOOST_CHECK(mesh.triangles().size() == 4);
      BOOST_CHECK(mesh.tetrahedra().size() == 1);
    }

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(1);
    Eigen::VectorXd readData(1);
    // Expected value in the middle of the tetra is the average of inputs (13.0/4)
    expected << 13.0 / 4;

    dt = participant.getMaxTimeStepSize();
    participant.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));

    participant.finalize();
  }
}

void testMappingVolumeOneTetraConservative(const std::string configFile, const TestContext &context)
{
  using precice::testing::equals;

  precice::Participant participant(context.name, context.config(), context.rank, context.size);
  // SolverOne defines a vertex and a conserved quantity (e.g. a force) on it.
  // SolverTwo defines a triangle and read the mapped quantity. We check it is spread correctly.

  std::vector<precice::VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords{0.1, 0.2, 0.3};
    vertexIDs.resize(coords.size() / 3);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    BOOST_TEST(vertexIDs[0] != -1, "Vertex A is invalid");
    BOOST_CHECK(participant.getMeshVertexSize(meshName) == 1);

    // Initialize, write data, advance and finalize
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values{1.0};
    participant.writeData(meshName, dataName, vertexIDs, values);

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant must advance only once.");
    participant.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords{0.0, 0.0, 0.0,
                               1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0};
    vertexIDs.resize(coords.size() / 3);

    participant.setMeshVertices(meshName, coords, vertexIDs);

    participant.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);

    BOOST_CHECK(participant.getMeshVertexSize(meshName) == 4);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshTwo");
    // setMeshTetrahedron currently adds underlying connectivity
    BOOST_REQUIRE(mesh.nVertices() == 4);
    BOOST_REQUIRE(mesh.edges().empty());
    BOOST_REQUIRE(mesh.triangles().empty());
    BOOST_REQUIRE(mesh.tetrahedra().size() == 1);

    BOOST_TEST(equals(mesh.tetrahedra()[0].getVolume(), 1.0 / 6), "Tetrahedron volume must be 1/6");

    // Initialize, read data, advance and finalize. Check expected mapping
    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    BOOST_REQUIRE(mesh.edges().size() == 6);
    BOOST_REQUIRE(mesh.triangles().size() == 4);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 1);
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant must advance once.");

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(4);
    Eigen::VectorXd readData(4);
    // For conservative load, each point takes a fraction of the load.
    // This fraction is the barycentric coordinate, and there the load is 1.
    // Input point is (0.1, 0.2, 0.3) and barycentric coordinates are thus (0.4, 0.1, 0.2, 0.3)
    expected << 0.4, 0.1, 0.2, 0.3;

    dt = participant.getMaxTimeStepSize();
    participant.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));

    participant.finalize();
  }
}

#endif

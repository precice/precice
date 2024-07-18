#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/precice.hpp"

void testMappingNearestProjection(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.3;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{1.0, 0.0, z};
  Vector3d coordOneC{1.0, 1.0, z};
  Vector3d coordOneD{0.0, 1.0, z};
  double   valOneA = 1.0;
  double   valOneB = 3.0;
  double   valOneC = 5.0;
  double   valOneD = 7.0;

  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z + 0.1};               // Maps to vertex A
  Vector3d coordTwoB{0.0, 0.5, z - 0.01};              // Maps to edge AD
  Vector3d coordTwoC{2.0 / 3.0, 1.0 / 3.0, z + 0.001}; // Maps to triangle ABC
  // This corresponds to the point C from mesh two on the triangle ABC on mesh one.
  Vector3d barycenterABC{0.3798734633239789, 0.24025307335204216, 0.3798734633239789};
  double   expectedValTwoA = 1.0;
  double   expectedValTwoB = 4.0;
  double   expectedValTwoC = Vector3d{valOneA, valOneB, valOneC}.dot(barycenterABC);

  if (context.isNamed("SolverOne")) {
    precice::Participant participant("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshOne";

    // Setup mesh one.
    int idA = participant.setMeshVertex(meshName, coordOneA);
    int idB = participant.setMeshVertex(meshName, coordOneB);
    int idC = participant.setMeshVertex(meshName, coordOneC);
    int idD = participant.setMeshVertex(meshName, coordOneD);

    if (defineEdgesExplicitly) {
      if (useBulkFunctions) {
        std::vector ids{idA, idB, idB, idC, idC, idD, idD, idA, idC, idA};
        participant.setMeshEdges(meshName, ids);
      } else {
        participant.setMeshEdge(meshName, idA, idB);
        participant.setMeshEdge(meshName, idB, idC);
        participant.setMeshEdge(meshName, idC, idD);
        participant.setMeshEdge(meshName, idD, idA);
        participant.setMeshEdge(meshName, idC, idA);
      }
    }

    if (useBulkFunctions) {
      std::vector ids{idA, idB, idC, idC, idD, idA};
      participant.setMeshTriangles(meshName, ids);
    } else {
      participant.setMeshTriangle(meshName, idA, idB, idC);
      participant.setMeshTriangle(meshName, idC, idD, idA);
    }

    // Initialize, thus sending the mesh.
    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    auto dataAID = "DataOne";
    BOOST_TEST(!participant.requiresGradientDataFor(meshName, dataAID));

    int    ids[]  = {idA, idB, idC, idD};
    double data[] = {valOneA, valOneB, valOneC, valOneD};
    participant.writeData(meshName, dataAID, ids, data);

    // Advance, thus send the data to the receiving partner.
    participant.advance(maxDt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant should have to advance once!");
    participant.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant participant("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshTwo";

    // Setup receiving mesh.
    int idA = participant.setMeshVertex(meshName, coordTwoA);
    int idB = participant.setMeshVertex(meshName, coordTwoB);
    int idC = participant.setMeshVertex(meshName, coordTwoC);

    // Initialize, thus receive the data and map.
    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto dataAID = "DataOne";
    BOOST_TEST(!participant.requiresGradientDataFor(meshName, dataAID));

    double value[3];
    int    ids[] = {idA, idB, idC};
    participant.readData(meshName, dataAID, ids, maxDt, value);

    BOOST_TEST(value[0] == expectedValTwoA);
    BOOST_TEST(value[1] == expectedValTwoB);
    BOOST_TEST(value[2] == expectedValTwoC);

    // Verify that there is only one time step necessary.
    participant.advance(maxDt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant should have to advance once!");
    participant.finalize();
  }
}

void testQuadMappingNearestProjection(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.3;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{1.0, 0.0, z};
  Vector3d coordOneC{0.999999999, 1.0, z}; // Forces diagonal 0-2 to be shorter.
  Vector3d coordOneD{0.0, 1.0, z};
  double   valOneA = 1.0;
  double   valOneB = 3.0;
  double   valOneC = 5.0;
  double   valOneD = 7.0;

  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z + 0.1};               // Maps to vertex A
  Vector3d coordTwoB{0.0, 0.5, z - 0.01};              // Maps to edge AD
  Vector3d coordTwoC{2.0 / 3.0, 1.0 / 3.0, z + 0.001}; // Maps to triangle ABC
  // This corresponds to the point C from mesh two on the triangle ABC on mesh one.
  Vector3d barycenterABC{0.3798734633239789, 0.24025307335204216, 0.3798734633239789};
  double   expectedValTwoA = 1.0;
  double   expectedValTwoB = 4.0;
  double   expectedValTwoC = Vector3d{valOneA, valOneB, valOneC}.dot(barycenterABC);

  if (context.isNamed("SolverOne")) {
    precice::Participant participant("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshOne";

    // Setup mesh one.
    int idA = participant.setMeshVertex(meshName, coordOneA);
    int idB = participant.setMeshVertex(meshName, coordOneB);
    int idC = participant.setMeshVertex(meshName, coordOneC);
    int idD = participant.setMeshVertex(meshName, coordOneD);

    if (defineEdgesExplicitly) {
      if (useBulkFunctions) {
        std::vector ids{idA, idB, idB, idC, idC, idD, idD, idA};
        participant.setMeshEdges(meshName, ids);
      } else {
        participant.setMeshEdge(meshName, idA, idB);
        participant.setMeshEdge(meshName, idB, idC);
        participant.setMeshEdge(meshName, idC, idD);
        participant.setMeshEdge(meshName, idD, idA);
      }
    }

    if (useBulkFunctions) {
      std::vector ids{idA, idB, idC, idD};
      participant.setMeshQuads(meshName, ids);
    } else {
      participant.setMeshQuad(meshName, idA, idB, idC, idD);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    BOOST_REQUIRE(mesh.nVertices() == 4);
    if (defineEdgesExplicitly) {
      BOOST_REQUIRE(mesh.edges().size() == 4);
    } else {
      BOOST_REQUIRE(mesh.edges().empty());
    }
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    // Initialize, thus sending the mesh.
    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();
    BOOST_TEST(mesh.edges().size() == 5);
    BOOST_TEST(mesh.triangles().size() == 2);

    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    auto   dataAID = "DataOne";
    int    ids[]   = {idA, idB, idC, idD};
    double data[]  = {valOneA, valOneB, valOneC, valOneD};
    participant.writeData(meshName, dataAID, ids, data);

    // Advance, thus send the data to the receiving partner.
    participant.advance(maxDt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant should have to advance once!");
    participant.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant participant("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshTwo";

    // Setup receiving mesh.
    int idA = participant.setMeshVertex(meshName, coordTwoA);
    int idB = participant.setMeshVertex(meshName, coordTwoB);
    int idC = participant.setMeshVertex(meshName, coordTwoC);

    // Initialize, thus receive the data and map.
    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto   dataAID = "DataOne";
    int    ids[]   = {idA, idB, idC};
    double values[3];
    participant.readData(meshName, dataAID, ids, maxDt, values);

    BOOST_TEST(values[0] == expectedValTwoA);
    BOOST_TEST(values[1] == expectedValTwoB);
    BOOST_TEST(values[2] == expectedValTwoC);

    // Verify that there is only one time step necessary.
    participant.advance(maxDt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant should have to advance once!");
    participant.finalize();
  }
}

void testQuadMappingNearestProjectionTallKite(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.0;

  // MeshOne
  Vector3d coordOneA{-0.2, 0.0, z};
  Vector3d coordOneB{0.0, 0.5, z};
  Vector3d coordOneC{0.2, 0.0, z};
  Vector3d coordOneD{0.0, -0.5, z};

  if (context.isNamed("SolverOne")) {
    precice::Participant participant("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshOne";

    // Setup mesh one.
    int idA = participant.setMeshVertex(meshName, coordOneA);
    int idB = participant.setMeshVertex(meshName, coordOneB);
    int idC = participant.setMeshVertex(meshName, coordOneC);
    int idD = participant.setMeshVertex(meshName, coordOneD);

    if (defineEdgesExplicitly) {
      if (useBulkFunctions) {
        std::vector ids{idA, idB, idB, idC, idC, idD, idD, idA};
        participant.setMeshEdges(meshName, ids);
      } else {
        participant.setMeshEdge(meshName, idA, idB);
        participant.setMeshEdge(meshName, idB, idC);
        participant.setMeshEdge(meshName, idC, idD);
        participant.setMeshEdge(meshName, idD, idA);
      }
    }

    if (useBulkFunctions) {
      std::vector ids{idA, idB, idC, idD};
      participant.setMeshQuads(meshName, ids);
    } else {
      participant.setMeshQuad(meshName, idA, idB, idC, idD);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    BOOST_REQUIRE(mesh.nVertices() == 4);
    if (defineEdgesExplicitly) {
      BOOST_REQUIRE(mesh.edges().size() == 4);
    } else {
      BOOST_REQUIRE(mesh.edges().empty());
    }
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    for (auto &edge : mesh.edges()) {
      BOOST_TEST(mesh::edgeLength(edge) < 0.6);
    }

    participant.finalize();
  }
}

void testQuadMappingNearestProjectionWideKite(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.0;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{0.5, 0.2, z};
  Vector3d coordOneC{1.0, 0.0, z};
  Vector3d coordOneD{0.5, -0.2, z};

  if (context.isNamed("SolverOne")) {
    Participant participant("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshName = "MeshOne";

    // Setup mesh one.
    int idA = participant.setMeshVertex(meshName, coordOneA);
    int idB = participant.setMeshVertex(meshName, coordOneB);
    int idC = participant.setMeshVertex(meshName, coordOneC);
    int idD = participant.setMeshVertex(meshName, coordOneD);

    if (defineEdgesExplicitly) {
      if (useBulkFunctions) {
        std::vector ids{idA, idB, idB, idC, idC, idD, idD, idA};
        participant.setMeshEdges(meshName, ids);
      } else {
        participant.setMeshEdge(meshName, idA, idB);
        participant.setMeshEdge(meshName, idB, idC);
        participant.setMeshEdge(meshName, idC, idD);
        participant.setMeshEdge(meshName, idD, idA);
      }
    }

    if (useBulkFunctions) {
      std::vector ids{idA, idB, idD, idC};
      participant.setMeshQuads(meshName, ids);
    } else {
      participant.setMeshQuad(meshName, idA, idB, idD, idC);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    BOOST_REQUIRE(mesh.nVertices() == 4);
    if (defineEdgesExplicitly) {
      BOOST_REQUIRE(mesh.edges().size() == 4);
    } else {
      BOOST_REQUIRE(mesh.edges().empty());
    }
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    for (auto &edge : mesh.edges()) {
      BOOST_TEST(mesh::edgeLength(edge) < 0.6);
    }

    participant.finalize();
  }
}

#endif

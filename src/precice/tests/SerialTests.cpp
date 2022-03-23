#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <fstream>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "action/RecorderAction.hpp"
#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "math/geometry.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/types.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

struct SerialTestFixture {

  std::string _pathToTests;

  SerialTestFixture()
  {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
  }
};

namespace {
std::vector<double> readDoublesFromTXTFile(const std::string &filename, int skip = 0)
{
  std::ifstream is{filename};
  if (skip > 0) {
    std::string ignore;
    while (skip--) {
      is >> ignore;
    }
  }
  return {std::istream_iterator<double>{is}, std::istream_iterator<double>{}};
}
} // namespace

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Serial, SerialTestFixture)

/// One solver uses incremental position set, read/write methods.
/// @todo This test uses resetmesh. How did this ever work?
#if 0
BOOST_AUTO_TEST_CASE(testExplicitWithDataExchange)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double counter = 0.0;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, _pathToTests + "explicit-mpi-single.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int meshOneID = cplInterface.getMeshID("MeshOne");
    /* int squareID = */ cplInterface.getMeshID("Test-Square");
    int              forcesID     = cplInterface.getDataID("Forces", meshOneID);
    int              velocitiesID = cplInterface.getDataID("Velocities", meshOneID);
    std::vector<int> indices(8);
    int              i = 0;

    //need one vertex to start
    Vector3d vertex = Vector3d::Zero();
    cplInterface.setMeshVertex(meshOneID, vertex.data());
    double maxDt = cplInterface.initialize();

    const auto &vertices = impl(cplInterface).mesh("Test-Square").vertices();
    while (cplInterface.isCouplingOngoing()) {
      impl(cplInterface).resetMesh(meshOneID);

      i = 0;
      for (auto &vertex : vertices) {
        int index = cplInterface.setMeshVertex(meshOneID, vertex.getCoords().data());
        BOOST_TEST(index == vertex.getID());
        indices.at(i) = index;
        i++;
      }
      //for (VertexIterator it = vertices.begin(); it != vertices.end(); it++) {
      for (auto &vertex : vertices) {
        Vector3d force(Vector3d::Constant(counter) + vertex.getCoords());
        cplInterface.writeVectorData(forcesID, vertex.getID(), force.data());
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()) {
        i = 0;
        for (auto &vertex : vertices) {
          Vector3d vel   = Vector3d::Zero();
          int      index = indices.at(i);
          i++;
          cplInterface.readVectorData(velocitiesID, index, vel.data());
          BOOST_TEST(vel == Vector3d::Constant(counter) + vertex.getCoords());
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    MeshID meshID = cplInterface.getMeshID("Test-Square");
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 1.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 1.0, 0.0).data());
    int    forcesID     = cplInterface.getDataID("Forces", meshID);
    int    velocitiesID = cplInterface.getDataID("Velocities", meshID);
    double maxDt        = cplInterface.initialize();
    auto & vertices     = impl(cplInterface).mesh("Test-Square").vertices();
    // SolverTwo does not start the coupled simulation and has, hence,
    // already received the first data to be validated.
    for (auto &vertex : vertices) {
      Vector3d force = Vector3d::Zero();
      cplInterface.readVectorData(forcesID, vertex.getID(), force.data());
      BOOST_TEST(force == Vector3d::Constant(counter) + vertex.getCoords());
    }
    counter += 1.0;

    while (cplInterface.isCouplingOngoing()) {
      for (auto &vertex : vertices) {
        Vector3d vel(Vector3d::Constant(counter - 1.0) + vertex.getCoords());
        cplInterface.writeVectorData(velocitiesID, vertex.getID(), vel.data());
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()) {
        for (auto &vertex : vertices) {
          Vector3d force = Vector3d::Zero();
          cplInterface.readVectorData(forcesID, vertex.getID(), force.data());
          BOOST_TEST(force == Vector3d::Constant(counter) + vertex.getCoords());
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}
#endif

/// One solver uses block set/get/read/write methods.
/// @todo This test uses resetmesh. How did this ever work?
#if 0
BOOST_AUTO_TEST_CASE(testExplicitWithBlockDataExchange)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double counter = 0.0;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, _pathToTests + "explicit-mpi-single-non-inc.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int             meshOneID      = cplInterface.getMeshID("MeshOne");
    double          maxDt          = cplInterface.initialize();
    int             forcesID       = cplInterface.getDataID("Forces", meshOneID);
    int             pressuresID    = cplInterface.getDataID("Pressures", meshOneID);
    int             velocitiesID   = cplInterface.getDataID("Velocities", meshOneID);
    int             temperaturesID = cplInterface.getDataID("Temperatures", meshOneID);
    auto &          vertices       = impl(cplInterface).mesh("Test-Square").vertices();
    int             size           = vertices.size();
    Eigen::VectorXd writePositions(size * 3);
    Eigen::VectorXd getWritePositions(size * 3);
    Eigen::VectorXd forces(size * 3);
    Eigen::VectorXd pressures(size);
    Eigen::VectorXi writeIDs(size);
    Eigen::VectorXi getWriteIDs(size);
    Eigen::VectorXd readPositions(size * 3);
    Eigen::VectorXd getReadPositions(size * 3);
    Eigen::VectorXd velocities(size * 3);
    Eigen::VectorXd temperatures(size);
    Eigen::VectorXd expectedVelocities(size * 3);
    Eigen::VectorXd expectedTemperatures(size);
    Eigen::VectorXi readIDs(size);
    Eigen::VectorXi getReadIDs(size);

    while (cplInterface.isCouplingOngoing()) {
      impl(cplInterface).resetMesh(meshOneID);
      for (auto &vertex : vertices) {
        for (int dim = 0; dim < 3; dim++) {
          writePositions(vertex.getID() * 3 + dim) = vertex.getCoords()(dim);
        }
      }
      cplInterface.setMeshVertices(meshOneID, size, writePositions.data(),
                                   writeIDs.data());
      for (auto &vertex : vertices) {
        // Vector3d force ( Vector3D(counter) + wrap<3,double>(vertex.getCoords()) );
        Vector3d force(Vector3d::Constant(counter) + vertex.getCoords());
        for (int dim = 0; dim < 3; dim++)
          forces(vertex.getID() * 3 + dim) = force(dim);
        pressures(vertex.getID()) = counter + vertex.getCoords()(0);
      }
      cplInterface.writeBlockVectorData(forcesID, size, writeIDs.data(), forces.data());
      cplInterface.writeBlockScalarData(pressuresID, size, writeIDs.data(), pressures.data());

      cplInterface.getMeshVertices(meshOneID, size, writeIDs.data(),
                                   getWritePositions.data());
      BOOST_TEST(writePositions == getWritePositions);

      cplInterface.getMeshVertexIDsFromPositions(meshOneID, size, writePositions.data(),
                                                 getWriteIDs.data());
      BOOST_TEST(writeIDs == getWriteIDs);
      //cplInterface.mapWrittenData(meshID);
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()) {
        for (auto &vertex : vertices) {
          for (int dim = 0; dim < 3; dim++) {
            int index                 = vertex.getID() * 3 + dim;
            readPositions(index)      = vertex.getCoords()(dim);
            expectedVelocities(index) = counter + vertex.getCoords()(dim);
          }
          expectedTemperatures(vertex.getID()) = counter + vertex.getCoords()(0);
        }
        impl(cplInterface).resetMesh(meshOneID);
        cplInterface.setMeshVertices(meshOneID, size, readPositions.data(), readIDs.data());
        cplInterface.mapReadDataTo(meshOneID);
        cplInterface.readBlockVectorData(velocitiesID, size, readIDs.data(),
                                         velocities.data());
        cplInterface.readBlockScalarData(temperaturesID, size, readIDs.data(),
                                         temperatures.data());
        BOOST_TEST(velocities == expectedVelocities);
        BOOST_TEST(temperatures == expectedTemperatures);

        counter += 1.0;
      }
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));

    int    squareID       = cplInterface.getMeshID("Test-Square");
    int    forcesID       = cplInterface.getDataID("Forces", squareID);
    int    pressuresID    = cplInterface.getDataID("Pressures", squareID);
    int    velocitiesID   = cplInterface.getDataID("Velocities", squareID);
    int    temperaturesID = cplInterface.getDataID("Temperatures", squareID);
    MeshID meshID         = cplInterface.getMeshID("Test-Square");
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 1.0, 0.0).data());
    cplInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 1.0, 0.0).data());
    double      maxDt    = cplInterface.initialize();
    const auto &vertices = impl(cplInterface).mesh("Test-Square").vertices();
    // SolverTwo does not start the coupled simulation and has, hence,
    // already received the first data to be validated.
    for (auto &vertex : vertices) {
      Vector3d force    = Vector3d::Zero();
      double   pressure = 0.0;
      cplInterface.readVectorData(forcesID, vertex.getID(), force.data());
      cplInterface.readScalarData(pressuresID, vertex.getID(), pressure);
      BOOST_TEST(force == Vector3d::Constant(counter) + vertex.getCoords());
      BOOST_TEST(pressure == counter + vertex.getCoords()(0));
    }
    counter += 1.0;

    while (cplInterface.isCouplingOngoing()) {
      for (auto &vertex : vertices) {
        Vector3d vel(Vector3d::Constant(counter - 1.0) + vertex.getCoords());
        cplInterface.writeVectorData(velocitiesID, vertex.getID(), vel.data());
        double temperature = counter - 1.0 + vertex.getCoords()(0);
        cplInterface.writeScalarData(temperaturesID, vertex.getID(), temperature);
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()) {
        for (auto &vertex : vertices) {
          Vector3d force    = Vector3d::Zero();
          double   pressure = 0.0;
          cplInterface.readVectorData(forcesID, vertex.getID(), force.data());
          cplInterface.readScalarData(pressuresID, vertex.getID(), pressure);
          BOOST_TEST(force == Vector3d::Constant(counter) + vertex.getCoords());
          BOOST_TEST(pressure == counter + vertex.getCoords()(0));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}
#endif




void testSummationAction(const std::string &configFile, TestContext const &context)
{
  using Eigen::Vector3d;

  if (context.isNamed("SolverTarget")) {
    // Expected values in the target solver
    double expectedValueA = 3.0;
    double expectedValueB = 7.0;
    double expectedValueC = 11.0;
    double expectedValueD = 15.0;

    // Target solver
    SolverInterface cplInterface(context.name, configFile, 0, 1);

    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    const MeshID meshID = cplInterface.getMeshID("MeshTarget");

    int idA = cplInterface.setMeshVertex(meshID, coordA.data());
    int idB = cplInterface.setMeshVertex(meshID, coordB.data());
    int idC = cplInterface.setMeshVertex(meshID, coordC.data());
    int idD = cplInterface.setMeshVertex(meshID, coordD.data());

    // Initialize, the mesh
    double dt = cplInterface.initialize();

    // Read the summed data from the mesh.
    int    dataAID = cplInterface.getDataID("Target", meshID);
    double valueA, valueB, valueC, valueD;

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.readScalarData(dataAID, idA, valueA);
      cplInterface.readScalarData(dataAID, idB, valueB);
      cplInterface.readScalarData(dataAID, idC, valueC);
      cplInterface.readScalarData(dataAID, idD, valueD);

      BOOST_TEST(valueA == expectedValueA);
      BOOST_TEST(valueB == expectedValueB);
      BOOST_TEST(valueC == expectedValueC);
      BOOST_TEST(valueD == expectedValueD);

      dt = cplInterface.advance(dt);
    }

    cplInterface.finalize();
  } else if (context.isNamed("SolverSourceOne")) {
    // Source solver one
    SolverInterface cplInterface(context.name, configFile, 0, 1);

    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    const MeshID meshID = cplInterface.getMeshID("MeshOne");

    int idA = cplInterface.setMeshVertex(meshID, coordA.data());
    int idB = cplInterface.setMeshVertex(meshID, coordB.data());
    int idC = cplInterface.setMeshVertex(meshID, coordC.data());
    int idD = cplInterface.setMeshVertex(meshID, coordD.data());

    // Initialize, the mesh
    double dt = cplInterface.initialize();

    int    dataAID = cplInterface.getDataID("SourceOne", meshID);
    double valueA  = 1.0;
    double valueB  = 3.0;
    double valueC  = 5.0;
    double valueD  = 7.0;

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataAID, idA, valueA);
      cplInterface.writeScalarData(dataAID, idB, valueB);
      cplInterface.writeScalarData(dataAID, idC, valueC);
      cplInterface.writeScalarData(dataAID, idD, valueD);

      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverSourceTwo"));
    // Source solver two
    SolverInterface cplInterface(context.name, configFile, 0, 1);
    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    const MeshID meshID = cplInterface.getMeshID("MeshTwo");

    int idA = cplInterface.setMeshVertex(meshID, coordA.data());
    int idB = cplInterface.setMeshVertex(meshID, coordB.data());
    int idC = cplInterface.setMeshVertex(meshID, coordC.data());
    int idD = cplInterface.setMeshVertex(meshID, coordD.data());

    // Initialize, the mesh
    double dt = cplInterface.initialize();

    int    dataAID = cplInterface.getDataID("SourceTwo", meshID);
    double valueA  = 2.0;
    double valueB  = 4.0;
    double valueC  = 6.0;
    double valueD  = 8.0;

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataAID, idA, valueA);
      cplInterface.writeScalarData(dataAID, idB, valueB);
      cplInterface.writeScalarData(dataAID, idC, valueC);
      cplInterface.writeScalarData(dataAID, idD, valueD);

      dt = cplInterface.advance(dt);
    }

    cplInterface.finalize();
  }
}

/**
 * @brief Test for summation action
 *
 */
BOOST_AUTO_TEST_CASE(testSummationActionTwoSources)
{
  PRECICE_TEST("SolverTarget"_on(1_rank), "SolverSourceOne"_on(1_rank), "SolverSourceTwo"_on(1_rank));
  const std::string configFile = _pathToTests + "summation-action.xml";
  testSummationAction(configFile, context);
}

void testWatchIntegral(const std::string &configFile, TestContext &context)
{
  using Eigen::Vector2d;

  if (context.isNamed("SolverOne")) {
    SolverInterface cplInterface(context.name, configFile, 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    const MeshID meshID = cplInterface.getMeshID("MeshOne");

    int idA = cplInterface.setMeshVertex(meshID, coordA.data());
    int idB = cplInterface.setMeshVertex(meshID, coordB.data());
    int idC = cplInterface.setMeshVertex(meshID, coordC.data());

    cplInterface.setMeshEdge(meshID, idA, idB);
    cplInterface.setMeshEdge(meshID, idB, idC);

    // Initialize, the mesh
    double dt = cplInterface.initialize();

    int    dataAID = cplInterface.getDataID("DataOne", meshID);
    double valueA  = 1.0;
    double valueB  = 2.0;
    double valueC  = 3.0;

    double increment = 1.0;

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataAID, idA, valueA);
      cplInterface.writeScalarData(dataAID, idB, valueB);
      cplInterface.writeScalarData(dataAID, idC, valueC);

      dt = cplInterface.advance(dt);

      valueA += increment;
      valueB += increment;
      valueC += increment;
    }
    cplInterface.finalize();
  } else if (context.isNamed("SolverTwo")) {

    SolverInterface cplInterface(context.name, configFile, 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    const int meshTwoID = cplInterface.getMeshID("MeshTwo");

    int idA = cplInterface.setMeshVertex(meshTwoID, coordA.data());
    int idB = cplInterface.setMeshVertex(meshTwoID, coordB.data());
    int idC = cplInterface.setMeshVertex(meshTwoID, coordC.data());

    cplInterface.setMeshEdge(meshTwoID, idA, idB);
    cplInterface.setMeshEdge(meshTwoID, idB, idC);

    // Initialize the mesh
    double dt = cplInterface.initialize();

    int    dataAID = cplInterface.getDataID("DataTwo", meshTwoID);
    double valueA, valueB, valueC;

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.readScalarData(dataAID, idA, valueA);
      cplInterface.readScalarData(dataAID, idB, valueB);
      cplInterface.readScalarData(dataAID, idC, valueC);

      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegral.log";
      auto        result   = readDoublesFromTXTFile(fileName, 4);
      auto        expected = std::vector<double>{
          1.0, 9.5, 0.0, 3.0,
          2.0, 12.5, 0.0, 3.0,
          3.0, 12.5, 0.0, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegralNoScale.log";
      auto        result   = readDoublesFromTXTFile(fileName, 4);
      auto        expected = std::vector<double>{
          1.0, 9.0, 0.0, 3.0,
          2.0, 12.0, 0.0, 3.0,
          3.0, 12.0, 0.0, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(testWatchIntegralScaleAndNoScale)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile = _pathToTests + "watch-integral.xml";
  testWatchIntegral(configFile, context);
}

void testQuadMappingNearestProjectionTallKite(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.0;

  // MeshOne
  Vector3d coordOneA{-0.2, 0.0, z};
  Vector3d coordOneB{0.0, 0.5, z};
  Vector3d coordOneC{0.2, 0.0, z};
  Vector3d coordOneD{0.0, -0.5, z};

  if (context.isNamed("SolverOne")) {
    SolverInterface cplInterface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = cplInterface.getMeshID("MeshOne");

    // Setup mesh one.
    int idA = cplInterface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = cplInterface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = cplInterface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = cplInterface.setMeshVertex(meshOneID, coordOneD.data());

    if (defineEdgesExplicitly) {

      int idAB = cplInterface.setMeshEdge(meshOneID, idA, idB);
      int idBC = cplInterface.setMeshEdge(meshOneID, idB, idC);
      int idCD = cplInterface.setMeshEdge(meshOneID, idC, idD);
      int idDA = cplInterface.setMeshEdge(meshOneID, idD, idA);

      cplInterface.setMeshQuad(meshOneID, idAB, idBC, idCD, idDA);

    } else {
      cplInterface.setMeshQuadWithEdges(meshOneID, idA, idB, idC, idD);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(cplInterface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 4);
    BOOST_REQUIRE(mesh.edges().size() == 5);
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    for (auto &edge : mesh.edges()) {
      BOOST_TEST(mesh::edgeLength(edge) < 0.6);
    }

    cplInterface.finalize();
  }
}

void testQuadMappingNearestProjectionWideKite(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.0;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{0.5, 0.2, z};
  Vector3d coordOneC{1.0, 0.0, z};
  Vector3d coordOneD{0.5, -0.2, z};

  if (context.isNamed("SolverOne")) {
    SolverInterface cplInterface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = cplInterface.getMeshID("MeshOne");

    // Setup mesh one.
    int idA = cplInterface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = cplInterface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = cplInterface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = cplInterface.setMeshVertex(meshOneID, coordOneD.data());

    if (defineEdgesExplicitly) {

      int idAB = cplInterface.setMeshEdge(meshOneID, idA, idB);
      int idBC = cplInterface.setMeshEdge(meshOneID, idB, idC);
      int idCD = cplInterface.setMeshEdge(meshOneID, idC, idD);
      int idDA = cplInterface.setMeshEdge(meshOneID, idD, idA);

      cplInterface.setMeshQuad(meshOneID, idAB, idCD, idBC, idDA);

    } else {
      cplInterface.setMeshQuadWithEdges(meshOneID, idA, idB, idD, idC);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(cplInterface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 4);
    BOOST_REQUIRE(mesh.edges().size() == 5);
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    for (auto &edge : mesh.edges()) {
      BOOST_TEST(mesh::edgeLength(edge) < 0.6);
    }

    cplInterface.finalize();
  }
}

void testQuadMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
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
    SolverInterface cplInterface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = cplInterface.getMeshID("MeshOne");

    // Setup mesh one.
    int idA = cplInterface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = cplInterface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = cplInterface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = cplInterface.setMeshVertex(meshOneID, coordOneD.data());

    if (defineEdgesExplicitly) {

      int idAB = cplInterface.setMeshEdge(meshOneID, idA, idB);
      int idBC = cplInterface.setMeshEdge(meshOneID, idB, idC);
      int idCD = cplInterface.setMeshEdge(meshOneID, idC, idD);
      int idDA = cplInterface.setMeshEdge(meshOneID, idD, idA);

      cplInterface.setMeshQuad(meshOneID, idAB, idBC, idCD, idDA);

    } else {
      cplInterface.setMeshQuadWithEdges(meshOneID, idA, idB, idC, idD);
    }

    auto &mesh = testing::WhiteboxAccessor::impl(cplInterface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 4);
    BOOST_REQUIRE(mesh.edges().size() == 5);
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    int dataAID = cplInterface.getDataID("DataOne", meshOneID);
    cplInterface.writeScalarData(dataAID, idA, valOneA);
    cplInterface.writeScalarData(dataAID, idB, valOneB);
    cplInterface.writeScalarData(dataAID, idC, valOneC);
    cplInterface.writeScalarData(dataAID, idD, valOneD);

    // Advance, thus send the data to the receiving partner.
    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    SolverInterface cplInterface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    int meshTwoID = cplInterface.getMeshID("MeshTwo");

    // Setup receiving mesh.
    int idA = cplInterface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = cplInterface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = cplInterface.setMeshVertex(meshTwoID, coordTwoC.data());

    // Initialize, thus receive the data and map.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    double valueA, valueB, valueC;
    cplInterface.readScalarData(dataAID, idA, valueA);
    cplInterface.readScalarData(dataAID, idB, valueB);
    cplInterface.readScalarData(dataAID, idC, valueC);

    BOOST_TEST(valueA == expectedValTwoA);
    BOOST_TEST(valueB == expectedValTwoB);
    BOOST_TEST(valueC == expectedValTwoC);

    // Verify that there is only one time step necessary.
    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    cplInterface.finalize();
  }
}

/**
 * @brief Tests the Nearest Projection Mapping between two participants with explicit definition of edges from a quad to a triangle
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingNearestProjectionExplicitEdges)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  bool              defineEdgesExplicitly = true;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjection(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping between two participants with explicit definition of edges from a quad to a triangle
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingNearestProjectionImplicitEdges)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  bool              defineEdgesExplicitly = false;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjection(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuad
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingDiagonalNearestProjectionExplicitEdgesTallKite)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  bool              defineEdgesExplicitly = true;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjectionTallKite(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuadWithEdges
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingDiagonalNearestProjectionImplicitEdgesTallKite)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  bool              defineEdgesExplicitly = false;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjectionTallKite(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuad
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingDiagonalNearestProjectionExplicitEdgesWideKite)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  bool              defineEdgesExplicitly = true;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjectionWideKite(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuadWithEdges
 *
 */
BOOST_AUTO_TEST_CASE(testQuadMappingDiagonalNearestProjectionImplicitEdgesWideKite)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  bool              defineEdgesExplicitly = false;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testQuadMappingNearestProjectionWideKite(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief method to test whether certain convergence measures give the correct number of iterations
 *
 */
void testConvergenceMeasures(const std::string configFile, TestContext const &context, std::vector<int> &expectedIterations)
{
  using Eigen::Vector2d;
  using namespace precice::constants;

  std::string meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";

  SolverInterface cplInterface(context.name, configFile, 0, 1);
  const int       meshID = cplInterface.getMeshID(meshName);

  Vector2d vertex{0.0, 0.0};

  std::vector<double> writeValues = {1.0, 1.01, 2.0, 2.5, 2.8, 2.81};

  VertexID vertexID = cplInterface.setMeshVertex(meshID, vertex.data());

  cplInterface.initialize();
  int numberOfAdvanceCalls = 0;
  int numberOfIterations   = -1;
  int timestep             = 0;

  while (cplInterface.isCouplingOngoing()) {
    if (cplInterface.isActionRequired(actionWriteIterationCheckpoint())) {
      numberOfIterations = 0;
      cplInterface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    if (context.isNamed("SolverTwo")) {
      DataID dataID = cplInterface.getDataID("Data2", meshID);
      cplInterface.writeScalarData(dataID, vertexID, writeValues.at(numberOfAdvanceCalls));
    }

    cplInterface.advance(1.0);
    ++numberOfAdvanceCalls;
    ++numberOfIterations;

    if (cplInterface.isActionRequired(actionReadIterationCheckpoint())) {
      cplInterface.markActionFulfilled(actionReadIterationCheckpoint());
    } else { //converged
      BOOST_TEST(numberOfIterations == expectedIterations.at(timestep));
      ++timestep;
    }
  }
  cplInterface.finalize();
}

BOOST_AUTO_TEST_CASE(testConvergenceMeasures1)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile         = _pathToTests + "convergence-measures1.xml";
  std::vector<int>  expectedIterations = {2, 4};
  testConvergenceMeasures(configFile, context, expectedIterations);
}

BOOST_AUTO_TEST_CASE(testConvergenceMeasures2)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile         = _pathToTests + "convergence-measures2.xml";
  std::vector<int>  expectedIterations = {3, 3};
  testConvergenceMeasures(configFile, context, expectedIterations);
}

BOOST_AUTO_TEST_CASE(testConvergenceMeasures3)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile         = _pathToTests + "convergence-measures3.xml";
  std::vector<int>  expectedIterations = {2, 4};
  testConvergenceMeasures(configFile, context, expectedIterations);
}



void testQuadMappingScaledConsistent(const std::string configFile, const TestContext &context)
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

  double expectedIntegral = math::geometry::triangleArea(coordOneA, coordOneB, coordOneC) * (valOneA + valOneB + valOneC) / 3.0 +
                            math::geometry::triangleArea(coordOneA, coordOneC, coordOneD) * (valOneA + valOneC + valOneD) / 3.0;

  if (context.isNamed("SolverOne")) {
    SolverInterface cplInterface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = cplInterface.getMeshID("MeshOne");

    // Setup mesh one.
    int idA = cplInterface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = cplInterface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = cplInterface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = cplInterface.setMeshVertex(meshOneID, coordOneD.data());

    int idAB = cplInterface.setMeshEdge(meshOneID, idA, idB);
    int idBC = cplInterface.setMeshEdge(meshOneID, idB, idC);
    int idCD = cplInterface.setMeshEdge(meshOneID, idC, idD);
    int idDA = cplInterface.setMeshEdge(meshOneID, idD, idA);

    cplInterface.setMeshQuad(meshOneID, idAB, idBC, idCD, idDA);

    auto &mesh = testing::WhiteboxAccessor::impl(cplInterface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 4);
    BOOST_REQUIRE(mesh.edges().size() == 5);
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    int dataAID = cplInterface.getDataID("DataOne", meshOneID);
    cplInterface.writeScalarData(dataAID, idA, valOneA);
    cplInterface.writeScalarData(dataAID, idB, valOneB);
    cplInterface.writeScalarData(dataAID, idC, valOneC);
    cplInterface.writeScalarData(dataAID, idD, valOneD);

    // Advance, thus send the data to the receiving partner.
    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    SolverInterface cplInterface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    int meshTwoID = cplInterface.getMeshID("MeshTwo");

    // Setup receiving mesh.
    int idA = cplInterface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = cplInterface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = cplInterface.setMeshVertex(meshTwoID, coordTwoC.data());

    int idAB = cplInterface.setMeshEdge(meshTwoID, idA, idB);
    int idBC = cplInterface.setMeshEdge(meshTwoID, idB, idC);
    int idAC = cplInterface.setMeshEdge(meshTwoID, idA, idC);

    cplInterface.setMeshTriangle(meshTwoID, idAB, idBC, idAC);

    // Initialize, thus receive the data and map.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    double valueA, valueB, valueC;
    cplInterface.readScalarData(dataAID, idA, valueA);
    cplInterface.readScalarData(dataAID, idB, valueB);
    cplInterface.readScalarData(dataAID, idC, valueC);

    double calculatedIntegral = math::geometry::triangleArea(coordTwoA, coordTwoB, coordTwoC) * (valueA + valueB + valueC) / 3.0;
    BOOST_TEST(expectedIntegral == calculatedIntegral);

    // Verify that there is only one time step necessary.
    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(testQuadMappingScaledConsistentOnA)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile = _pathToTests + "mapping-scaled-consistent-onA.xml";
  testQuadMappingScaledConsistent(configFile, context);
}

BOOST_AUTO_TEST_CASE(testQuadMappingScaledConsistentOnB)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile = _pathToTests + "mapping-scaled-consistent-onB.xml";
  testQuadMappingScaledConsistent(configFile, context);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

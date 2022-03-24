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


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

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
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

struct SerialTestFixture : testing::WhiteboxAccessor {

  std::string _pathToTests;

  void reset()
  {
    mesh::Data::resetDataCount();
  }

  SerialTestFixture()
  {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
    reset();
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

/// Test reading of a full features coupling configuration file.
BOOST_AUTO_TEST_CASE(TestConfigurationPeano)
{
  PRECICE_TEST(1_rank);
  std::string filename = _pathToTests + "configuration.xml";

  // Test configuration for accessor "Peano"
  SolverInterface interfacePeano("Peano", filename, 0, 1);

  BOOST_TEST(impl(interfacePeano)._participants.size() == 2);
  BOOST_TEST(interfacePeano.getDimensions() == 2);

  impl::PtrParticipant peano = impl(interfacePeano)._participants.at(0);
  BOOST_TEST(peano);
  BOOST_TEST(peano->getName() == "Peano");

  std::vector<impl::MeshContext *> meshContexts = peano->_meshContexts;
  BOOST_TEST(meshContexts.size() == 2);
  BOOST_TEST(peano->_usedMeshContexts.size() == 2);

  BOOST_TEST(meshContexts.at(0)->mesh->getName() == std::string("PeanoNodes"));
  BOOST_TEST(meshContexts.at(1)->mesh->getName() == std::string("ComsolNodes"));
}

BOOST_AUTO_TEST_CASE(TestConfigurationComsol)
{
  PRECICE_TEST(1_rank);
  std::string filename = _pathToTests + "configuration.xml";

  // Test configuration for accessor "Comsol"
  SolverInterface interfaceComsol("Comsol", filename, 0, 1);
  BOOST_TEST(impl(interfaceComsol)._participants.size() == 2);
  BOOST_TEST(interfaceComsol.getDimensions() == 2);

  impl::PtrParticipant comsol = impl(interfaceComsol)._participants.at(1);
  BOOST_TEST(comsol);
  BOOST_TEST(comsol->getName() == "Comsol");

  std::vector<impl::MeshContext *> meshContexts = comsol->_meshContexts;
  BOOST_TEST(meshContexts.size() == 2);
  BOOST_TEST(meshContexts.at(0) == static_cast<void *>(nullptr));
  BOOST_TEST(meshContexts.at(1)->mesh->getName() == std::string("ComsolNodes"));
  BOOST_TEST(comsol->_usedMeshContexts.size() == 1);
}

BOOST_AUTO_TEST_SUITE(Lifecycle)

// Test representing the full explicit lifecycle of a SolverInterface
BOOST_AUTO_TEST_CASE(Full)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshid   = interface.getMeshID("MeshOne");
    double coords[] = {0.1, 1.2, 2.3};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto   dataid = interface.getDataID("DataOne", meshid);
    double data[] = {3.4, 4.5, 5.6};
    interface.writeVectorData(dataid, vertexid, data);
  } else {
    auto   meshid   = interface.getMeshID("MeshTwo");
    double coords[] = {0.12, 1.21, 2.2};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto dataid = interface.getDataID("DataTwo", meshid);
    interface.writeScalarData(dataid, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
  interface.finalize();
}

// Test representing the full lifecycle of a SolverInterface
// Finalize is not called explicitly here.
// The destructor has to cleanup.
BOOST_AUTO_TEST_CASE(ImplicitFinalize)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshid   = interface.getMeshID("MeshOne");
    double coords[] = {0.1, 1.2, 2.3};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto   dataid = interface.getDataID("DataOne", meshid);
    double data[] = {3.4, 4.5, 5.6};
    interface.writeVectorData(dataid, vertexid, data);
  } else {
    auto   meshid   = interface.getMeshID("MeshTwo");
    double coords[] = {0.12, 1.21, 2.2};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto dataid = interface.getDataID("DataTwo", meshid);
    interface.writeScalarData(dataid, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
}

// Test representing the minimal lifecylce, which consists out of construction only.
// The destructor has to cleanup correctly.
BOOST_AUTO_TEST_CASE(ConstructOnly)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);
}

// Test representing the minimal lifecylce with explicit finalization.
// This shows how to manually finalize MPI etc without using the SolverInterface.
BOOST_AUTO_TEST_CASE(ConstructAndExplicitFinalize)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END()

/// Test to run simple "do nothing" coupling between two solvers.
void runTestExplicit(std::string const &configurationFileName, TestContext const &context)
{
  BOOST_TEST_MESSAGE("Config: " << configurationFileName);

  int    timesteps = 0;
  double time      = 0.0;

  SolverInterface couplingInterface(context.name, configurationFileName, 0, 1);

  //was necessary to replace pre-defined geometries
  if (context.isNamed("SolverOne") && couplingInterface.hasMesh("MeshOne")) {
    int meshID = couplingInterface.getMeshID("MeshOne");
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
  }
  if (context.isNamed("SolverTwo") && couplingInterface.hasMesh("Test-Square")) {
    int meshID = couplingInterface.getMeshID("Test-Square");
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
  }

  BOOST_TEST(couplingInterface.getDimensions() == 3);
  double dt = couplingInterface.initialize();
  while (couplingInterface.isCouplingOngoing()) {
    time += dt;
    dt = couplingInterface.advance(dt);
    timesteps++;
  }
  couplingInterface.finalize();

  BOOST_TEST(time == 10.0);
  BOOST_TEST(timesteps == 10);
}

BOOST_AUTO_TEST_CASE(TestExplicitMPISingle)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "explicit-mpi-single.xml";
  runTestExplicit(config, context);
}

BOOST_AUTO_TEST_CASE(TestExplicitMPI)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "explicit-mpi.xml";
  runTestExplicit(config, context);
}

BOOST_AUTO_TEST_CASE(TestExplicitSockets)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::string config = _pathToTests + "explicit-sockets.xml";
  runTestExplicit(config, context);
}

/// Test to run a simple "do nothing" coupling with subcycling solvers.
BOOST_AUTO_TEST_CASE(testExplicitWithSubcycling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, _pathToTests + "explicit-mpi-single.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    double maxDt     = precice.initialize();
    int    timestep  = 0;
    double dt        = maxDt / 2.0; // Timestep length desired by solver
    double currentDt = dt;          // Timestep length used by solver
    while (precice.isCouplingOngoing()) {
      maxDt     = precice.advance(currentDt);
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    BOOST_TEST(timestep == 20);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int meshID = precice.getMeshID("Test-Square");
    precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    precice.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    double maxDt     = precice.initialize();
    int    timestep  = 0;
    double dt        = maxDt / 3.0; // Timestep length desired by solver
    double currentDt = dt;          // Timestep length used by solver
    while (precice.isCouplingOngoing()) {
      maxDt     = precice.advance(currentDt);
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    BOOST_TEST(timestep == 30);
  }
}

/// One solver uses incremental position set, read/write methods.
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
    int meshID = cplInterface.getMeshID("Test-Square");
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

/**
 * @brief The second solver initializes the data of the first.
 *
 * A mapping is employed for the second solver, i.e., at the end of
 * initializeData(), the mapping needs to be invoked.
 */
BOOST_AUTO_TEST_CASE(testExplicitWithDataInitialization)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, _pathToTests + "explicit-data-init.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int meshOneID = cplInterface.getMeshID("MeshOne");
    cplInterface.setMeshVertex(meshOneID, Vector3d(1.0, 2.0, 3.0).data());
    double maxDt      = cplInterface.initialize();
    int    dataAID    = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID    = cplInterface.getDataID("DataTwo", meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(2.5 == valueDataB);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshTwoID, pos.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);
    cplInterface.writeScalarData(dataBID, 0, 2.0);
    //sagen dass daten jetzt geschrieben
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 2.5);
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

/// One solver uses block set/get/read/write methods.
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

    int squareID       = cplInterface.getMeshID("Test-Square");
    int forcesID       = cplInterface.getDataID("Forces", squareID);
    int pressuresID    = cplInterface.getDataID("Pressures", squareID);
    int velocitiesID   = cplInterface.getDataID("Velocities", squareID);
    int temperaturesID = cplInterface.getDataID("Temperatures", squareID);
    int meshID         = cplInterface.getMeshID("Test-Square");
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

/**
  * @brief Runs a coupled simulation where one solver supplies a geometry.
  *
  * SolverOne only reads the displacements of the geometry and checks whether
  * they are equals to the coordinates of SolverTwo. SolverTwo creates and
  * displaces the coordinates.
  *
  * @todo Maybe remove this test.
  */
BOOST_AUTO_TEST_CASE(testExplicitWithSolverGeometry)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  int    timesteps = 0;
  double time      = 0;

  SolverInterface couplingInterface(context.name, _pathToTests + "explicit-solvergeometry.xml", 0, 1);
  BOOST_TEST(couplingInterface.getDimensions() == 3);
  if (context.isNamed("SolverOne")) {
    //was necessary to replace pre-defined geometries
    int meshID = couplingInterface.getMeshID("MeshOne");
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());

    double dt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int meshID = couplingInterface.getMeshID("SolverGeometry");
    int i0     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    int i1     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    int i2     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 1.0, 0.0).data());
    int e0     = couplingInterface.setMeshEdge(meshID, i0, i1);
    int e1     = couplingInterface.setMeshEdge(meshID, i1, i2);
    int e2     = couplingInterface.setMeshEdge(meshID, i2, i0);
    couplingInterface.setMeshTriangle(meshID, e0, e1, e2);
    double dt = couplingInterface.initialize();

    int size = couplingInterface.getMeshVertexSize(meshID);
    BOOST_TEST(size == 3);

    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();
  }
}

/**
 * @brief Runs a coupled sim. with data scaling applied.
 *
 * SolverOne writes vector data on a cube geometry. The data values are defined
 * and stay constant over the coupling cycles. SolverTwo has a scaling of the
 * values activated and reads the scaled values.
 */
BOOST_AUTO_TEST_CASE(testExplicitWithDataScaling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface cplInterface(context.name, _pathToTests + "explicit-datascaling.xml", 0, 1);
  BOOST_TEST(cplInterface.getDimensions() == 2);

  std::vector<double> positions = {0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.0};
  std::vector<int>    ids       = {0, 0, 0, 0};

  if (context.isNamed("SolverOne")) {
    int meshID = cplInterface.getMeshID("Test-Square-One");
    cplInterface.setMeshVertices(meshID, 4, positions.data(), ids.data());
    for (int i = 0; i < 4; i++)
      cplInterface.setMeshEdge(meshID, ids.at(i), ids.at((i + 1) % 4));

    double dt = cplInterface.initialize();

    int velocitiesID = cplInterface.getDataID("Velocities", meshID);
    while (cplInterface.isCouplingOngoing()) {
      for (size_t i = 0; i < impl(cplInterface).mesh("Test-Square-One").vertices().size(); ++i) {
        Eigen::Vector2d data = Eigen::Vector2d::Constant(i);
        cplInterface.writeVectorData(velocitiesID, i, data.data());
      }
      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int meshID = cplInterface.getMeshID("Test-Square-Two");
    cplInterface.setMeshVertices(meshID, 4, positions.data(), ids.data());
    for (int i = 0; i < 4; i++)
      cplInterface.setMeshEdge(meshID, ids.at(i), ids.at((i + 1) % 4));

    double dt = cplInterface.initialize();

    int velocitiesID = cplInterface.getDataID("Velocities", meshID);
    while (cplInterface.isCouplingOngoing()) {
      const auto size = impl(cplInterface).mesh("Test-Square-Two").vertices().size();
      for (size_t i = 0; i < size; ++i) {
        Eigen::Vector2d readData;
        cplInterface.readVectorData(velocitiesID, i, readData.data());
        Eigen::Vector2d expectedData = Eigen::Vector2d::Constant(i * 10.0);
        BOOST_TEST(readData(0) == expectedData(0));
        BOOST_TEST(readData(1) == expectedData(1));
      }
      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();
  }
}

/// Test simple coupled simulation with coupling iterations.
BOOST_AUTO_TEST_CASE(testImplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double state              = 0.0;
  double checkpoint         = 0.0;
  int    iterationCount     = 0;
  double initialStateChange = 5.0;
  double stateChange        = initialStateChange;
  int    computedTimesteps  = 0;
  using namespace precice::constants;

  SolverInterface couplingInterface(context.name, _pathToTests + "implicit.xml", 0, 1);

  if (context.isNamed("SolverOne")) {
    int    meshID = couplingInterface.getMeshID("Square");
    double pos[3];
    // Set mesh positions
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    couplingInterface.setMeshVertex(meshID, pos);
    pos[0] = 1.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    couplingInterface.setMeshVertex(meshID, pos);
    pos[0] = 1.0;
    pos[1] = 1.0;
    pos[2] = 0.0;
    couplingInterface.setMeshVertex(meshID, pos);
    pos[0] = 0.0;
    pos[1] = 1.0;
    pos[2] = 0.0;
    couplingInterface.setMeshVertex(meshID, pos);

    double maxDt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())) {
        couplingInterface.markActionFulfilled(actionWriteIterationCheckpoint());
        checkpoint     = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())) {
        couplingInterface.markActionFulfilled(actionReadIterationCheckpoint());
        state = checkpoint;
      }
      iterationCount++;
      stateChange = initialStateChange / (double) iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimeWindowComplete()) {
        computedTimesteps++;
      }
    }
    couplingInterface.finalize();
    BOOST_TEST(computedTimesteps == 4);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    double maxDt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())) {
        couplingInterface.markActionFulfilled(actionWriteIterationCheckpoint());
        checkpoint     = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())) {
        couplingInterface.markActionFulfilled(actionReadIterationCheckpoint());
        state = checkpoint;
        iterationCount++;
      }
      stateChange = initialStateChange / (double) iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimeWindowComplete()) {
        computedTimesteps++;
      }
    }
    couplingInterface.finalize();
    BOOST_TEST(computedTimesteps == 4);
  }
}

/// Test simple coupled simulation with iterations, data initialization and without acceleration
BOOST_AUTO_TEST_CASE(testImplicitWithDataInitialization)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using namespace precice::constants;

  SolverInterface couplingInterface(context.name, _pathToTests + "implicit-data-init.xml", 0, 1);

  int         dimensions = couplingInterface.getDimensions();
  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  double      writeValue, expectedReadValue;

  if (context.isNamed("SolverOne")) {
    meshName          = "MeshOne";
    writeDataName     = "Forces";
    readDataName      = "Velocities";
    writeValue        = 1;
    expectedReadValue = 2;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName          = "MeshTwo";
    writeDataName     = "Velocities";
    readDataName      = "Forces";
    writeValue        = 2;
    expectedReadValue = 1;
  }
  int                 meshID      = couplingInterface.getMeshID(meshName);
  int                 writeDataID = couplingInterface.getDataID(writeDataName, meshID);
  int                 readDataID  = couplingInterface.getDataID(readDataName, meshID);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshID, vertex.data());

  double dt = 0;
  dt        = couplingInterface.initialize();
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);
  const std::string & cowid = actionWriteInitialData();

  if (couplingInterface.isActionRequired(cowid)) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    couplingInterface.markActionFulfilled(cowid);
  }

  couplingInterface.initializeData();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())) {
      couplingInterface.markActionFulfilled(actionWriteIterationCheckpoint());
    }
    couplingInterface.readVectorData(readDataID, vertexID, readData.data());
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    dt = couplingInterface.advance(dt);
    if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())) {
      couplingInterface.markActionFulfilled(actionReadIterationCheckpoint());
    }
  }
  couplingInterface.finalize();
}

/// Tests stationary mapping with solver provided meshes.
void runTestStationaryMappingWithSolverMesh(std::string const &config, int dim, TestContext const &context)
{
  std::string meshForcesA = "MeshForcesA";
  std::string meshDisplA  = "MeshDisplacementsA";
  std::string meshForcesB = "MeshForcesB";
  std::string meshDisplB  = "MeshDisplacementsB";
  std::string dataForces  = "Forces";
  std::string dataDispl   = "Displacements";
  using testing::equals;

  SolverInterface interface(context.name, config, 0, 1);
  BOOST_TEST(interface.getDimensions() == dim);

  std::vector<Eigen::VectorXd> positions;
  Eigen::VectorXd              position(dim);
  if (dim == 2) {
    position << 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 0.0;
    positions.push_back(position);
    position << 1.0, 1.0;
    positions.push_back(position);
    position << 0.0, 1.0;
    positions.push_back(position);
  } else {
    position << 0.0, 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 1.0, 0.0;
    positions.push_back(position);
    position << 0.0, 1.0, 1.0;
    positions.push_back(position);
    position << 0.0, 0.0, 1.0;
    positions.push_back(position);
  }
  size_t size = positions.size();

  if (context.isNamed("SolverA")) {
    int meshForcesID = interface.getMeshID(meshForcesA);
    int meshDisplID  = interface.getMeshID(meshDisplA);
    int dataForcesID = interface.getDataID(dataForces, meshForcesID);
    int dataDisplID  = interface.getDataID(dataDispl, meshDisplID);

    // Set solver mesh positions for reading and writing data with mappings
    for (size_t i = 0; i < size; i++) {
      position = positions.at(i).array() + 0.1;
      interface.setMeshVertex(meshForcesID, position.data());
      position = positions.at(i).array() + 0.6;
      interface.setMeshVertex(meshDisplID, position.data());
    }
    double maxDt = interface.initialize();

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(not interface.isReadDataAvailable());
    Eigen::VectorXd force = Eigen::VectorXd::Constant(dim, 1);
    Eigen::VectorXd displ = Eigen::VectorXd::Constant(dim, 0);
    for (size_t i = 0; i < size; i++) {
      interface.writeVectorData(dataForcesID, i, force.data());
    }
    interface.mapWriteDataFrom(meshForcesID);
    maxDt = interface.advance(maxDt);
    interface.mapReadDataTo(meshDisplID);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    force.array() += 1.0;
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataDisplID, i, displ.data());
      BOOST_TEST(displ(0) == positions.at(i)(0) + 0.1);
      interface.writeVectorData(dataForcesID, i, force.data());
    }
    interface.mapWriteDataFrom(meshForcesID);
    maxDt = interface.advance(maxDt);
    interface.mapReadDataTo(meshDisplID);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataDisplID, i, displ.data());
      BOOST_TEST(displ(0) == 2.0 * (positions.at(i)(0) + 0.1));
    }
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverB"));
    int meshForcesID = interface.getMeshID(meshForcesB);
    int meshDisplID  = interface.getMeshID(meshDisplB);
    int dataForcesID = interface.getDataID(dataForces, meshForcesID);
    int dataDisplID  = interface.getDataID(dataDispl, meshDisplID);

    // Set solver mesh positions provided to SolverA for data mapping
    for (size_t i = 0; i < size; i++) {
      interface.setMeshVertex(meshForcesID, positions.at(i).data());
      position = positions.at(i).array() + 0.5;
      interface.setMeshVertex(meshDisplID, position.data());
    }
    double maxDt = interface.initialize();

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    Eigen::VectorXd force      = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd totalForce = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd displ      = Eigen::VectorXd::Zero(dim);
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataForcesID, i, force.data());
      totalForce += force;
      displ.setConstant(positions.at(i)(0));
      interface.writeVectorData(dataDisplID, i, displ.data());
    }
    Eigen::VectorXd expected = Eigen::VectorXd::Constant(dim, size);
    BOOST_TEST(equals(totalForce, expected));
    maxDt = interface.advance(maxDt);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    totalForce.setConstant(0);
    for (size_t i = 0; i < positions.size(); i++) {
      interface.readVectorData(dataForcesID, i, force.data());
      totalForce += force;
      displ.setConstant(2.0 * positions.at(i)(0));
      interface.writeVectorData(dataDisplID, i, displ.data());
    }
    expected.setConstant(2.0 * (double) size);
    BOOST_TEST(equals(totalForce, expected));
    maxDt = interface.advance(maxDt);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(not interface.isReadDataAvailable()); //second participant has no new data after last advance
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataForcesID, i, force.data());
    }
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(testStationaryMappingWithSolverMesh2D)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank));
  std::string config = _pathToTests + "mapping-without-geo-2D.xml";
  runTestStationaryMappingWithSolverMesh(config, 2, context);
}

BOOST_AUTO_TEST_CASE(testStationaryMappingWithSolverMesh3D)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank));
  std::string config = _pathToTests + "mapping-without-geo-3D.xml";
  runTestStationaryMappingWithSolverMesh(config, 3, context);
}

/**
 * @brief Buggy simulation setup of FSI coupling between Flite and Calculix.
 *
 * Bug: after first call of advance by Flite the mapped forces are value NaN.
 *
 * Some information about the coupling:
 * - explicit coupling scheme
 * - Flite (incompressible Navier-Stokes) starts simulation
 * - Mapping is done on Flite side with RBF
 *
 * @todo rename this test and config
 */
BOOST_AUTO_TEST_CASE(testBug)
{
  PRECICE_TEST("Flite"_on(1_rank), "Calculix"_on(1_rank));

  using Eigen::Vector3d;
  std::string configName = _pathToTests + "bug.xml";

  int                   slices = 5;
  std::vector<Vector3d> coords;
  for (int i = 0; i < slices; i++) {
    double z = (double) i * 1.0;
    coords.push_back(Vector3d(1.0, 0.0, z));
    coords.push_back(Vector3d(0.0, 1.0, z));
    coords.push_back(Vector3d(-1.0, 0.0, z));
    coords.push_back(Vector3d(0.0, -1.0, z));
  }

  if (context.isNamed("Flite")) {
    SolverInterface precice("Flite", configName, 0, 1);

    int meshID             = precice.getMeshID("FliteNodes");
    int forcesID           = precice.getDataID("Forces", meshID);
    int displacementsID    = precice.getDataID("Displacements", meshID);
    int oldDisplacementsID = precice.getDataID("OldDisplacements", meshID);
    BOOST_TEST(precice.getDimensions() == 3);
    for (Vector3d &coord : coords) {
      precice.setMeshVertex(meshID, coord.data());
    }
    double maxDt = precice.initialize();
    double dt    = 1.0e-5 / 15.0; // Flite took 15 subcycling steps
    while (precice.isCouplingOngoing()) {
      dt = dt < maxDt ? dt : maxDt;
      for (int i = 0; i < (int) coords.size(); i++) {
        double force[3] = {1.0, 2.0, 3.0};
        precice.writeVectorData(forcesID, i, force);
      }
      maxDt = precice.advance(dt);
      precice.mapReadDataTo(meshID);
      for (int i = 0; i < (int) coords.size(); i++) {
        double displacement[3];
        double oldDisplacement[3];
        precice.readVectorData(displacementsID, i, displacement);
        precice.readVectorData(oldDisplacementsID, i, oldDisplacement);
      }
    }
    precice.finalize();
  } else {
    BOOST_TEST(context.isNamed("Calculix"));
    SolverInterface precice("Calculix", configName, 0, 1);

    int meshID = precice.getMeshID("CalculixNodes");
    for (Vector3d &coord : coords) {
      precice.setMeshVertex(meshID, coord.data());
    }
    for (int i = 0; i < slices - 1; i++) {
      // Build cylinder/channel geometry
      precice.setMeshTriangleWithEdges(meshID, i * 4, (i * 4) + 1, (i + 1) * 4);
      precice.setMeshTriangleWithEdges(meshID, (i + 1) * 4, (i * 4) + 1, ((i + 1) * 4) + 1);
      precice.setMeshTriangleWithEdges(meshID, i * 4 + 1, (i * 4) + 2, (i + 1) * 4 + 1);
      precice.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 1, (i * 4) + 2, ((i + 1) * 4) + 2);
      precice.setMeshTriangleWithEdges(meshID, i * 4 + 2, (i * 4) + 3, (i + 1) * 4 + 2);
      precice.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 2, (i * 4) + 3, ((i + 1) * 4) + 3);
      precice.setMeshTriangleWithEdges(meshID, i * 4 + 3, (i * 4), (i + 1) * 4 + 3);
      precice.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 3, i * 4, (i + 1) * 4);
    }
    double dt = precice.initialize();
    while (precice.isCouplingOngoing()) {
      precice.advance(dt);
    }
    precice.finalize();
  }
}

/**
 * @brief Three solvers are coupled in a fork S2 <-> S1 <-> S3.
 *
 * Both couplings are explicit, solver 1 provides the mesh to the other two
 * solvers.
 */
void runTestThreeSolvers(std::string const &config, std::vector<int> expectedCallsOfAdvance, TestContext const &context)
{
  std::string writeIterCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterCheckpoint(constants::actionReadIterationCheckpoint());
  std::string writeInitData(constants::actionWriteInitialData());

  int callsOfAdvance = 0;

  if (context.isNamed("SolverOne")) {
    SolverInterface precice(context.name, config, 0, 1);

    int meshAID = precice.getMeshID("MeshA");
    int meshBID = precice.getMeshID("MeshB");
    precice.setMeshVertex(meshAID, Eigen::Vector2d(0, 0).data());
    precice.setMeshVertex(meshBID, Eigen::Vector2d(1, 1).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(0));
  } else if (context.isNamed("SolverTwo")) {
    SolverInterface precice(context.name, config, 0, 1);

    int meshID = precice.getMeshID("MeshC");
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(1));
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    SolverInterface precice(context.name, config, 0, 1);

    int meshID = precice.getMeshID("MeshD");
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(2));
  }
}

BOOST_AUTO_TEST_CASE(ThreeSolversExplicitExplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string      config = _pathToTests + "three-solver-explicit-explicit.xml";
  std::vector<int> expectedCallsOfAdvance{10, 10, 10};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

BOOST_AUTO_TEST_CASE(ThreeSolversImplicitImplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string      config = _pathToTests + "three-solver-implicit-implicit.xml";
  std::vector<int> expectedCallsOfAdvance{30, 30, 20};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

BOOST_AUTO_TEST_CASE(ThreeSolversImplicitExplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string      config = _pathToTests + "three-solver-implicit-explicit.xml";
  std::vector<int> expectedCallsOfAdvance{30, 30, 10};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

BOOST_AUTO_TEST_CASE(ThreeSolversExplicitImplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string      config = _pathToTests + "three-solver-explicit-implicit.xml";
  std::vector<int> expectedCallsOfAdvance{30, 10, 30};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

BOOST_AUTO_TEST_CASE(ThreeSolversParallel)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string      config = _pathToTests + "three-solver-parallel.xml";
  std::vector<int> expectedCallsOfAdvance{30, 30, 10};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

/// Four solvers are multi-coupled.
BOOST_AUTO_TEST_CASE(MultiCoupling)
{
  PRECICE_TEST("SOLIDZ1"_on(1_rank), "SOLIDZ2"_on(1_rank), "SOLIDZ3"_on(1_rank), "NASTIN"_on(1_rank));
  std::vector<Eigen::Vector2d> positions;
  Eigen::Vector2d              position;
  position << 0.0, 0.0;
  positions.push_back(position);
  position << 1.0, 0.0;
  positions.push_back(position);
  position << 1.0, 1.0;
  positions.push_back(position);
  position << 0.0, 1.0;
  positions.push_back(position);

  std::vector<Eigen::Vector2d> datas;
  Eigen::Vector2d              data;
  data << 1.0, 1.0;
  datas.push_back(data);
  data << 2.0, 2.0;
  datas.push_back(position);
  data << 3.0, 3.0;
  datas.push_back(data);
  data << 4.0, 5.0;
  datas.push_back(data);

  std::string writeIterCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterCheckpoint(constants::actionReadIterationCheckpoint());

  if (context.isNamed("SOLIDZ1") ||
      context.isNamed("SOLIDZ2") ||
      context.isNamed("SOLIDZ3")) {
    int meshID      = -1;
    int dataWriteID = -1;
    int dataReadID  = -1;

    SolverInterface precice(context.name, _pathToTests + "/multi.xml", 0, 1);
    BOOST_TEST(precice.getDimensions() == 2);

    if (context.isNamed("SOLIDZ1")) {
      meshID      = precice.getMeshID("SOLIDZ_Mesh1");
      dataWriteID = precice.getDataID("Displacements1", meshID);
      dataReadID  = precice.getDataID("Forces1", meshID);
    } else if (context.isNamed("SOLIDZ2")) {
      meshID      = precice.getMeshID("SOLIDZ_Mesh2");
      dataWriteID = precice.getDataID("Displacements2", meshID);
      dataReadID  = precice.getDataID("Forces2", meshID);
    } else if (context.isNamed("SOLIDZ3")) {
      meshID      = precice.getMeshID("SOLIDZ_Mesh3");
      dataWriteID = precice.getDataID("Displacements3", meshID);
      dataReadID  = precice.getDataID("Forces3", meshID);
    }

    std::vector<int> vertexIDs;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshID, positions.at(i).data());
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i = 0; i < 4; i++) {
      precice.writeVectorData(dataWriteID, vertexIDs.at(i), datas.at(i).data());
    }

    if (precice.isActionRequired(writeIterCheckpoint)) {
      precice.markActionFulfilled(writeIterCheckpoint);
    }
    precice.advance(0.0001);
    if (precice.isActionRequired(readIterCheckpoint)) {
      precice.markActionFulfilled(readIterCheckpoint);
    }

    for (size_t i = 0; i < 4; i++) {
      precice.readVectorData(dataReadID, vertexIDs.at(i), datas.at(i).data());
    }

    BOOST_TEST(datas.at(0)(0) == 1.00000000000000002082e-03);
    BOOST_TEST(datas.at(0)(1) == 1.00000000000000002082e-03);
    BOOST_TEST(datas.at(1)(0) == 0.00000000000000000000e+00);
    BOOST_TEST(datas.at(1)(1) == 1.00000000000000002082e-03);
    BOOST_TEST(datas.at(2)(0) == 3.00000000000000006245e-03);
    BOOST_TEST(datas.at(2)(1) == 3.00000000000000006245e-03);
    BOOST_TEST(datas.at(3)(0) == 4.00000000000000008327e-03);
    BOOST_TEST(datas.at(3)(1) == 5.00000000000000010408e-03);

    precice.finalize();

  } else {
    BOOST_TEST(context.isNamed("NASTIN"));
    SolverInterface precice("NASTIN", _pathToTests + "/multi.xml", 0, 1);
    BOOST_TEST(precice.getDimensions() == 2);
    int meshID1      = precice.getMeshID("NASTIN_Mesh1");
    int meshID2      = precice.getMeshID("NASTIN_Mesh2");
    int meshID3      = precice.getMeshID("NASTIN_Mesh3");
    int dataWriteID1 = precice.getDataID("Forces1", meshID1);
    int dataWriteID2 = precice.getDataID("Forces2", meshID2);
    int dataWriteID3 = precice.getDataID("Forces3", meshID3);

    std::vector<int> vertexIDs1;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshID1, positions.at(i).data());
      vertexIDs1.push_back(vertexID);
    }
    std::vector<int> vertexIDs2;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshID2, positions.at(i).data());
      vertexIDs2.push_back(vertexID);
    }
    std::vector<int> vertexIDs3;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshID3, positions.at(i).data());
      vertexIDs3.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i = 0; i < 4; i++) {
      precice.writeVectorData(dataWriteID1, vertexIDs1.at(i), datas.at(i).data());
      precice.writeVectorData(dataWriteID2, vertexIDs2.at(i), datas.at(i).data());
      precice.writeVectorData(dataWriteID3, vertexIDs3.at(i), datas.at(i).data());
    }

    if (precice.isActionRequired(writeIterCheckpoint)) {
      precice.markActionFulfilled(writeIterCheckpoint);
    }
    precice.advance(0.0001);
    if (precice.isActionRequired(readIterCheckpoint)) {
      precice.markActionFulfilled(readIterCheckpoint);
    }

    precice.finalize();
  }
}

void testMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
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
      int idCA = cplInterface.setMeshEdge(meshOneID, idC, idA);

      cplInterface.setMeshTriangle(meshOneID, idAB, idBC, idCA);
      cplInterface.setMeshTriangle(meshOneID, idCD, idDA, idCA);

    } else {
      cplInterface.setMeshTriangleWithEdges(meshOneID, idA, idB, idC);
      cplInterface.setMeshTriangleWithEdges(meshOneID, idC, idD, idA);
    }

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
 * @brief Tests the Nearest Projection Mapping between two participants with explicit definition of edges
 *
 */
BOOST_AUTO_TEST_CASE(MappingNearestProjectionExplicitEdges)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  bool              defineEdgesExplicitly = true;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testMappingNearestProjection(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests the Nearest Projection Mapping between two participants with explicit definition of edges
 *
 */
BOOST_AUTO_TEST_CASE(MappingNearestProjectionImplicitEdges)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  bool              defineEdgesExplicitly = false;
  const std::string configFile            = _pathToTests + "mapping-nearest-projection.xml";
  testMappingNearestProjection(defineEdgesExplicitly, configFile, context);
}

/**
 * @brief Tests sending one mesh to multiple participants
 *
 */
BOOST_AUTO_TEST_CASE(SendMeshToMultipleParticipants)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));

  const std::string configFile = _pathToTests + "send-mesh-to-multiple-participants.xml";
  std::string       meshName;

  Eigen::Vector2d vertex{0.0, 0.0};

  double value = 1.0;

  if (context.isNamed("SolverOne")) {
    meshName = "MeshA";
  } else if (context.isNamed("SolverTwo")) {
    meshName = "MeshB";
  } else if (context.isNamed("SolverThree")) {
    meshName = "MeshC";
  }

  SolverInterface cplInterface(context.name, configFile, 0, 1);

  const int meshID = cplInterface.getMeshID(meshName);

  int vertexID = cplInterface.setMeshVertex(meshID, vertex.data());

  double maxDt = cplInterface.initialize();

  int dataID = cplInterface.getDataID("Data", meshID);

  if (context.isNamed("SolverOne")) {
    cplInterface.writeScalarData(dataID, vertexID, value);
  } else {
    double valueReceived = -1.0;
    cplInterface.readScalarData(dataID, vertexID, valueReceived);
    BOOST_TEST(valueReceived == value);
  }

  cplInterface.advance(maxDt);
  cplInterface.finalize();
}

/**
 * @brief Test to reproduce the problem of issue 383, https://github.com/precice/precice/issues/383
 *
 */
BOOST_AUTO_TEST_CASE(PreconditionerBug)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  const std::string configFile = _pathToTests + "preconditioner-bug.xml";

  std::string meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";

  SolverInterface cplInterface(context.name, configFile, 0, 1);
  const int       meshID = cplInterface.getMeshID(meshName);

  Vector2d vertex{0.0, 0.0};

  int vertexID = cplInterface.setMeshVertex(meshID, vertex.data());

  cplInterface.initialize();
  int numberOfAdvanceCalls = 0;

  while (cplInterface.isCouplingOngoing()) {
    if (cplInterface.isActionRequired(actionWriteIterationCheckpoint()))
      cplInterface.markActionFulfilled(actionWriteIterationCheckpoint());
    if (cplInterface.isActionRequired(actionReadIterationCheckpoint()))
      cplInterface.markActionFulfilled(actionReadIterationCheckpoint());

    if (context.isNamed("SolverTwo")) {
      int dataID = cplInterface.getDataID("DataOne", meshID);
      // to get convergence in first timestep (everything 0), but not in second timestep
      Vector2d value{0.0, 2.0 + numberOfAdvanceCalls * numberOfAdvanceCalls};
      cplInterface.writeVectorData(dataID, vertexID, value.data());
    }
    cplInterface.advance(1.0);
    ++numberOfAdvanceCalls;
  }
  cplInterface.finalize();
}

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

    const int meshID = cplInterface.getMeshID("MeshTarget");

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

    const int meshID = cplInterface.getMeshID("MeshOne");

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

    const int meshID = cplInterface.getMeshID("MeshTwo");

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

    const int meshID = cplInterface.getMeshID("MeshOne");

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

  int vertexID = cplInterface.setMeshVertex(meshID, vertex.data());

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
      int dataID = cplInterface.getDataID("Data2", meshID);
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

/**
 * @brief Test to make sure that actions are called in the right order for explicit coupling via RecorderAction
 */
BOOST_AUTO_TEST_CASE(ActionTimingsExplicit)
{
  const std::string configFile = _pathToTests + "recorder-action-explicit.xml";
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using namespace precice::constants;

  SolverInterface couplingInterface(context.name, configFile, 0, 1);

  int         dimensions = couplingInterface.getDimensions();
  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  double      writeValue;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Forces";
    readDataName  = "Velocities";
    writeValue    = 1;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Velocities";
    readDataName  = "Forces";
    writeValue    = 2;
  }
  int                 meshID      = couplingInterface.getMeshID(meshName);
  int                 writeDataID = couplingInterface.getDataID(writeDataName, meshID);
  int                 readDataID  = couplingInterface.getDataID(readDataName, meshID);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshID, vertex.data());

  double dt = -1;
  BOOST_TEST(action::RecorderAction::records.empty());
  dt = couplingInterface.initialize();
  BOOST_TEST(dt == 1.0);
  if (context.isNamed("SolverOne")) {
    BOOST_TEST(action::RecorderAction::records.empty());
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    BOOST_TEST(action::RecorderAction::records.size() == 2);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::READ_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::READ_MAPPING_POST);
  }
  action::RecorderAction::reset();
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);
  const std::string & cowid = actionWriteInitialData();

  if (couplingInterface.isActionRequired(cowid)) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    couplingInterface.markActionFulfilled(cowid);
  }

  couplingInterface.initializeData();
  if (context.isNamed("SolverOne")) {
    BOOST_TEST(action::RecorderAction::records.size() == 2);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    BOOST_TEST(action::RecorderAction::records.size() == 4);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
    BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::READ_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(3).timing == action::Action::READ_MAPPING_POST);
  }
  action::RecorderAction::reset();

  int iteration = 0;

  while (couplingInterface.isCouplingOngoing()) {
    couplingInterface.readVectorData(readDataID, vertexID, readData.data());
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    dt = couplingInterface.advance(dt);
    BOOST_TEST(couplingInterface.isTimeWindowComplete());
    iteration++;
    if (context.isNamed("SolverOne") || iteration < 10) {
      BOOST_TEST(action::RecorderAction::records.size() == 5);
      BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
      BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::READ_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(3).timing == action::Action::READ_MAPPING_POST);
      BOOST_TEST(action::RecorderAction::records.at(4).timing == action::Action::ON_TIME_WINDOW_COMPLETE_POST);
    } else { // SolverTwo only writes in very last iteration, does not read.
      BOOST_TEST(action::RecorderAction::records.size() == 3);
      BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
      BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::ON_TIME_WINDOW_COMPLETE_POST);
    }
    action::RecorderAction::reset();
  }
  couplingInterface.finalize();
}

/**
 * @brief Test to make sure that actions are called in the right order for implicit coupling via RecorderAction
 */
BOOST_AUTO_TEST_CASE(ActionTimingsImplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  const std::string configFile = _pathToTests + "recorder-action-implicit.xml";

  using namespace precice::constants;

  SolverInterface couplingInterface(context.name, configFile, 0, 1);

  int         dimensions = couplingInterface.getDimensions();
  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  std::string writeIterCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterCheckpoint(constants::actionReadIterationCheckpoint());
  double      writeValue;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Forces";
    readDataName  = "Velocities";
    writeValue    = 1;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Velocities";
    readDataName  = "Forces";
    writeValue    = 2;
  }
  int                 meshID      = couplingInterface.getMeshID(meshName);
  int                 writeDataID = couplingInterface.getDataID(writeDataName, meshID);
  int                 readDataID  = couplingInterface.getDataID(readDataName, meshID);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshID, vertex.data());

  double dt = -1;
  BOOST_TEST(action::RecorderAction::records.empty());
  dt = couplingInterface.initialize();
  BOOST_TEST(dt == 1.0);
  if (context.isNamed("SolverOne")) {
    BOOST_TEST(action::RecorderAction::records.empty());
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    BOOST_TEST(action::RecorderAction::records.size() == 2);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::READ_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::READ_MAPPING_POST);
  }
  action::RecorderAction::reset();
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);
  const std::string & cowid = actionWriteInitialData();

  if (couplingInterface.isActionRequired(cowid)) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    couplingInterface.markActionFulfilled(cowid);
  }

  couplingInterface.initializeData();
  if (context.isNamed("SolverOne")) {
    BOOST_TEST(action::RecorderAction::records.size() == 2);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    BOOST_TEST(action::RecorderAction::records.size() == 4);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
    BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::READ_MAPPING_PRIOR);
    BOOST_TEST(action::RecorderAction::records.at(3).timing == action::Action::READ_MAPPING_POST);
  }
  action::RecorderAction::reset();

  int iteration = 0;

  while (couplingInterface.isCouplingOngoing()) {
    couplingInterface.readVectorData(readDataID, vertexID, readData.data());
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    if (couplingInterface.isActionRequired(writeIterCheckpoint)) {
      couplingInterface.markActionFulfilled(writeIterCheckpoint);
    }
    dt = couplingInterface.advance(dt);
    if (couplingInterface.isActionRequired(readIterCheckpoint)) {
      couplingInterface.markActionFulfilled(readIterCheckpoint);
    }
    if (couplingInterface.isTimeWindowComplete()) {
      iteration++;
    }
    if (context.isNamed("SolverOne") || iteration < 10) {
      if (couplingInterface.isTimeWindowComplete()) {
        BOOST_TEST(action::RecorderAction::records.size() == 5);
        BOOST_TEST(action::RecorderAction::records.at(4).timing == action::Action::ON_TIME_WINDOW_COMPLETE_POST);
      } else {
        BOOST_TEST(action::RecorderAction::records.size() == 4);
      }
      BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
      BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::READ_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(3).timing == action::Action::READ_MAPPING_POST);
    } else { // SolverTwo only writes in very last iteration, does not read.
      if (couplingInterface.isTimeWindowComplete()) {
        BOOST_TEST(action::RecorderAction::records.size() == 3);
        BOOST_TEST(action::RecorderAction::records.at(2).timing == action::Action::ON_TIME_WINDOW_COMPLETE_POST);
      } else {
        BOOST_TEST(action::RecorderAction::records.size() == 2);
      }
      BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_PRIOR);
      BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::WRITE_MAPPING_POST);
    }
    action::RecorderAction::reset();
  }
  couplingInterface.finalize();
}

BOOST_AUTO_TEST_CASE(MultipleFromMappings)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  const std::string configFile = _pathToTests + "multiple-from-mappings.xml";

  SolverInterface interface(context.name, configFile, 0, 1);
  Vector2d        vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const int meshIDTop      = interface.getMeshID("MeshATop");
    const int meshIDBottom   = interface.getMeshID("MeshABottom");
    int       vertexIDTop    = interface.setMeshVertex(meshIDTop, vertex.data());
    int       vertexIDBottom = interface.setMeshVertex(meshIDBottom, vertex.data());
    int       dataIDTop      = interface.getDataID("Pressure", meshIDTop);
    int       dataIDBottom   = interface.getDataID("Pressure", meshIDBottom);

    double dt = interface.initialize();
    interface.advance(dt);
    double pressure = -1.0;
    interface.readScalarData(dataIDTop, vertexIDTop, pressure);
    BOOST_TEST(pressure == 1.0);
    pressure = -1.0;
    interface.readScalarData(dataIDBottom, vertexIDBottom, pressure);
    BOOST_TEST(pressure == 1.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const int meshID   = interface.getMeshID("MeshB");
    int       vertexID = interface.setMeshVertex(meshID, vertex.data());
    int       dataID   = interface.getDataID("Pressure", meshID);

    double dt       = interface.initialize();
    double pressure = 1.0;
    interface.writeScalarData(dataID, vertexID, pressure);
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(MultipleToMappings)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  const std::string configFile = _pathToTests + "multiple-to-mappings.xml";

  SolverInterface interface(context.name, configFile, 0, 1);
  Vector2d        vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const int meshIDTop      = interface.getMeshID("MeshATop");
    const int meshIDBottom   = interface.getMeshID("MeshABottom");
    int       vertexIDTop    = interface.setMeshVertex(meshIDTop, vertex.data());
    int       vertexIDBottom = interface.setMeshVertex(meshIDBottom, vertex.data());
    int       dataIDTop      = interface.getDataID("DisplacementTop", meshIDTop);
    int       dataIDBottom   = interface.getDataID("DisplacementBottom", meshIDBottom);

    double dt              = interface.initialize();
    double displacementTop = 1.0;
    interface.writeScalarData(dataIDTop, vertexIDTop, displacementTop);
    double displacementBottom = 2.0;
    interface.writeScalarData(dataIDBottom, vertexIDBottom, displacementBottom);
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const int meshID   = interface.getMeshID("MeshB");
    int       vertexID = interface.setMeshVertex(meshID, vertex.data());
    int       dataID   = interface.getDataID("DisplacementSum", meshID);

    double dt = interface.initialize();
    interface.advance(dt);
    double displacement = -1.0;
    interface.readScalarData(dataID, vertexID, displacement);
    BOOST_TEST(displacement == 3.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(AitkenAcceleration)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  const std::string configFile = _pathToTests + "aitken-acceleration.xml";

  SolverInterface interface(context.name, configFile, 0, 1);
  Vector2d        vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const int meshID   = interface.getMeshID("A-Mesh");
    int       vertexID = interface.setMeshVertex(meshID, vertex.data());
    int       dataID   = interface.getDataID("Data", meshID);

    double dt    = interface.initialize();
    double value = 1.0;
    interface.writeScalarData(dataID, vertexID, value);

    interface.markActionFulfilled(actionWriteIterationCheckpoint());
    interface.advance(dt);
    interface.markActionFulfilled(actionReadIterationCheckpoint());
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const int meshID   = interface.getMeshID("B-Mesh");
    int       vertexID = interface.setMeshVertex(meshID, vertex.data());
    int       dataID   = interface.getDataID("Data", meshID);

    double dt = interface.initialize();
    interface.markActionFulfilled(actionWriteIterationCheckpoint());
    interface.advance(dt);

    double value = -1.0;
    interface.readScalarData(dataID, vertexID, value);
    BOOST_TEST(value == 0.1); // due to initial underrelaxation

    interface.markActionFulfilled(actionReadIterationCheckpoint());
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

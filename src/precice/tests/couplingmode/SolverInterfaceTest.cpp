#include "SolverInterfaceTest.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/config/Configuration.hpp"
#include "geometry/Geometry.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include <fstream>

#include "tarch/tests/TestCaseFactory.h"
registerIntegrationTest(precice::tests::SolverInterfaceTest)

namespace precice {
namespace tests {

logging::Logger SolverInterfaceTest:: _log("tests::SolverInterfaceTest");

SolverInterfaceTest:: SolverInterfaceTest()
:
  tarch::tests::TestCase("tests::SolverInterfaceTest"),
  _pathToTests()
{}

void SolverInterfaceTest:: setUp()
{
  _pathToTests = utils::Globals::getPathToSources() + "/precice/tests/couplingmode/";
}

void SolverInterfaceTest:: run()
{
  TRACE();
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testConfiguration();
  }
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 1){
    MPI_Comm comm = Par::getRestrictedCommunicator({0, 1});
    if (Par::getProcessRank() <= 1){
      Par::setGlobalCommunicator(comm);
      testMethod(testExplicit);
      testMethod(testExplicitWithSubcycling);
      testMethod(testExplicitWithDataExchange);
      testMethod(testExplicitWithDataInitialization);
      testMethod(testExplicitWithBlockDataExchange);
      testMethod(testExplicitWithSolverGeometry);
      testMethod(testExplicitWithDisplacingGeometry);
      testMethod(testExplicitWithDataScaling);
#     ifndef PRECICE_NO_SPIRIT2
      testMethod(testExplicitWithCheckpointingStatMapping);
#     endif // not PRECICE_NO_SPIRIT2
      testMethod(testImplicit);
#     ifndef PRECICE_NO_SPIRIT2
      testMethod(testImplicitWithCheckpointingMappingStat);
#     endif // not PRECICE_NO_SPIRIT2
      testMethod(testStationaryMappingWithSolverMesh);
      testMethod(testBug);
      testMethod(testNASTINMeshRestart);
#     ifndef PRECICE_NO_PYTHON
      testMethod(testPinelliCoupled);
#     endif
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  if (Par::getCommunicatorSize() > 2){
    MPI_Comm comm = Par::getRestrictedCommunicator({0, 1, 2});
    if (Par::getProcessRank() <= 2){
      Par::setGlobalCommunicator(comm);
      testMethod(testThreeSolvers);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  Par::synchronizeProcesses();
  if (Par::getCommunicatorSize() > 3){
    MPI_Comm comm = Par::getRestrictedCommunicator({0, 1, 2, 3});
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod(testDistributedCommunications)
      testMethod(testMultiCoupling);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

#ifndef PRECICE_NO_MPI

void SolverInterfaceTest:: configureSolverInterface
(
  const std::string& configFilename,
  SolverInterface&   interface )
{
  TRACE(configFilename);
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  precice::utils::EventRegistry::clear();
  utils::MasterSlave::reset();
  config::Configuration config;
  utils::configure(config.getXMLTag(), configFilename);
  //validate ( config.isValid() );
  interface._impl->configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceTest:: testConfiguration()
{
  TRACE();
  std::string filename = _pathToTests + "/configuration.xml";
  // Test configuration for accessor "Peano"
  SolverInterface interfacePeano ("Peano", 0, 1);
  configureSolverInterface ( filename, interfacePeano );
  validateEquals(interfacePeano._impl->_participants.size(), 2 );
  validateEquals(interfacePeano.getDimensions(), 2 );

  impl::PtrParticipant peano = interfacePeano._impl->_participants[0];
  validate(peano.use_count() > 0 );
  validateEquals(peano->getName(), "Peano");
  validateEquals(peano->getID(), 0);

  std::vector<impl::MeshContext*> meshContexts = peano->_meshContexts;
  validateEquals ( meshContexts.size(), 3 );
  validateEquals ( peano->_usedMeshContexts.size(), 3 );

  validateEquals ( meshContexts[0]->mesh->getName(), std::string("PeanoNodes") );
  validateEquals ( meshContexts[1]->mesh->getName(), std::string("ComsolNodes") );
  validateEquals ( meshContexts[2]->mesh->getName(), std::string("Boundingbox") );

  // Test configuration for accessor "Comsol"
  SolverInterface interfaceComsol ("Comsol", 0, 1);
  configureSolverInterface ( filename, interfaceComsol );
  validateEquals ( interfaceComsol._impl->_participants.size(), 2 );
  validateEquals ( interfaceComsol.getDimensions(), 2 );

  impl::PtrParticipant comsol = interfaceComsol._impl->_participants[1];
  validate ( comsol.use_count() > 0 );
  validateEquals ( comsol->getName(), "Comsol" );
  validateEquals ( comsol->getID(), 1 );

  meshContexts = comsol->_meshContexts;
  validateEquals(meshContexts.size(), 3);
  // Can be replaced by nullptr as soon as we have C++11 available.
  validateEquals(meshContexts[0], static_cast<void*>(nullptr));
  validateEquals(meshContexts[1]->mesh->getName(), std::string("ComsolNodes"));
    // Can be replaced by nullptr as soon as we have C++11 available.
  validateEquals(meshContexts[2], static_cast<void*>(nullptr));
  validateEquals(comsol->_usedMeshContexts.size(), 1);
}

void SolverInterfaceTest:: testExplicit()
{
  TRACE();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );

  int timesteps;
  double time;

  if (utils::Parallel::getProcessRank() == 0){
    runSolver("SolverOne", "/explicit-mpi-single.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    runSolver("SolverTwo", "/explicit-mpi-single.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }

  if ( utils::Parallel::getProcessRank() == 0){
    runSolver("SolverOne", "/explicit-mpi.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    runSolver("SolverTwo", "/explicit-mpi.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }

  #ifndef PRECICE_NO_SOCKETS
  if (utils::Parallel::getProcessRank() == 0){
    runSolver("SolverOne", "/explicit-sockets.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    runSolver("SolverTwo", "/explicit-sockets.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
  utils::Parallel::synchronizeProcesses(); // close all sockets before continuing
  #endif // not PRECICE_NO_SOCKETS

  if (utils::Parallel::getProcessRank() == 0){
    runSolver("SolverOne", "/explicit-files.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    runSolver("SolverTwo", "/explicit-files.xml",
              timesteps, time);
    validateEquals(time, 10.0);
    validateEquals(timesteps, 10);
  }
}

void SolverInterfaceTest:: testExplicitWithSubcycling()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);

  mesh::Mesh::resetGeometryIDsGlobally();

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface precice("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-mpi-single.xml", precice);
    double maxDt = precice.initialize();
    int timestep = 0;
    double dt = maxDt / 2.0; // Timestep length desired by solver
    double currentDt = dt;   // Timestep length used by solver
    while (precice.isCouplingOngoing()){
      maxDt = precice.advance(currentDt);
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    validateEquals(timestep, 20);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface precice("SolverTwo", 0, 1);
    configureSolverInterface ( _pathToTests + "/explicit-mpi-single.xml", precice );
    double maxDt = precice.initialize ();
    int timestep = 0;
    double dt = maxDt / 3.0; // Timestep length desired by solver
    double currentDt = dt;   // Timestep length used by solver
    while ( precice.isCouplingOngoing() ){
      maxDt = precice.advance ( currentDt );
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    validateEquals ( timestep, 30 );
  }
}

void SolverInterfaceTest:: testExplicitWithDataExchange()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  double counter = 0.0;
  using Eigen::Vector3d;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface cplInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-mpi-single.xml", cplInterface);

    int meshOneID = cplInterface.getMeshID("MeshOne");
    /* int squareID = */ cplInterface.getMeshID("Test-Square");
    int forcesID = cplInterface.getDataID("Forces", meshOneID);
    int velocitiesID = cplInterface.getDataID("Velocities", meshOneID);
    int indices[8];
    int i = 0;

    //need one vertex to start
    Vector3d vertex = Vector3d::Zero();
    cplInterface.setMeshVertex(meshOneID, vertex.data());
    double maxDt = cplInterface.initialize();

    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    while (cplInterface.isCouplingOngoing()){
      cplInterface.resetMesh(meshOneID);
      i = 0;
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        int index = cplInterface.setMeshVertex(meshOneID, it.vertexCoords());
        validateEquals(index, it.vertexID());
        indices[i] = index;
        i++;
      }
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3d force(Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords()));
        cplInterface.writeVectorData(forcesID, it.vertexID(), force.data());
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        i=0;
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          Vector3d vel = Vector3d::Zero();
          int index = indices[i];
          i++;
          cplInterface.readVectorData(velocitiesID, index, vel.data());
          validate(math::equals(vel, Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords())));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface cplInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-mpi-single.xml", cplInterface);
    double maxDt = cplInterface.initialize();
    int meshID = cplInterface.getMeshID("Test-Square");
    int forcesID = cplInterface.getDataID("Forces", meshID);
    int velocitiesID = cplInterface.getDataID("Velocities", meshID);
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    // SolverTwo does not start the coupled simulation and has, hence,
    // already received the first data to be validated.
    for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
      Vector3d force = Vector3d::Zero();
      cplInterface.readVectorData(forcesID, it.vertexID(), force.data());
      validate(math::equals(force, Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords())));
    }
    counter += 1.0;

    while (cplInterface.isCouplingOngoing()){
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3d vel(Vector3d::Constant(counter - 1.0) + Eigen::Map<const Vector3d>(it.vertexCoords()));
        cplInterface.writeVectorData(velocitiesID, it.vertexID(), vel.data());
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          Vector3d force = Vector3d::Zero();
          cplInterface.readVectorData(forcesID, it.vertexID(), force.data());
          validate(math::equals(force, Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords())));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithDataInitialization()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  using Eigen::Vector3d;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface cplInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-data-init.xml", cplInterface);
    int meshOneID = cplInterface.getMeshID("MeshOne");
    cplInterface.setMeshVertex(meshOneID, Vector3d(1.0,2.0,3.0).data());
    double maxDt = cplInterface.initialize();
    int dataAID = cplInterface.getDataID("DataOne",meshOneID);
    int dataBID = cplInterface.getDataID("DataTwo",meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    validateNumericalEquals(2.0, valueDataB);
    while (cplInterface.isCouplingOngoing()){
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readScalarData(dataBID, 0, valueDataB);
      validateNumericalEquals(2.5, valueDataB);
    }
    cplInterface.finalize();
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface cplInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-data-init.xml", cplInterface);
    int meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos = Vector3d::Zero();
    cplInterface.setMeshVertex(meshTwoID, pos.data());
    double maxDt = cplInterface.initialize();
    int dataAID = cplInterface.getDataID("DataOne",meshTwoID);
    int dataBID = cplInterface.getDataID("DataTwo",meshTwoID);
    cplInterface.writeScalarData(dataBID, 0, 2.0);
    //sagen dass daten jetzt geschrieben
    cplInterface.fulfilledAction(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    validateWithParams2(math::equals(valueDataA, expected), valueDataA, expected);
    while (cplInterface.isCouplingOngoing()){
      cplInterface.writeScalarData(dataBID, 0, 2.5);
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      validateWithParams2(math::equals(valueDataA, expected), valueDataA, expected);
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithBlockDataExchange()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  double counter = 0.0;
  using Eigen::Vector3d;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface cplInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-mpi-single-non-inc.xml",
                             cplInterface);
    double maxDt = cplInterface.initialize();
    int meshOneID = cplInterface.getMeshID("MeshOne");
    int forcesID = cplInterface.getDataID("Forces", meshOneID);
    int pressuresID = cplInterface.getDataID("Pressures", meshOneID);
    int velocitiesID = cplInterface.getDataID("Velocities", meshOneID);
    int temperaturesID = cplInterface.getDataID("Temperatures", meshOneID);
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    int size = vertices.size();
    Eigen::VectorXd writePositions(size*3);
    Eigen::VectorXd getWritePositions(size*3);
    Eigen::VectorXd forces(size*3);
    Eigen::VectorXd pressures(size);
    Eigen::VectorXi writeIDs(size);
    Eigen::VectorXi getWriteIDs(size);
    Eigen::VectorXd readPositions(size*3);
    Eigen::VectorXd getReadPositions(size*3);
    Eigen::VectorXd velocities(size*3);
    Eigen::VectorXd temperatures(size);
    Eigen::VectorXd expectedVelocities(size*3);
    Eigen::VectorXd expectedTemperatures(size);
    Eigen::VectorXi readIDs(size);
    Eigen::VectorXi getReadIDs(size);

    while (cplInterface.isCouplingOngoing()){
      cplInterface.resetMesh(meshOneID);
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        for (int dim=0; dim < 3; dim++){
          writePositions[it.vertexID()*3 + dim] = it.vertexCoords()[dim];
        }
      }
      cplInterface.setMeshVertices(meshOneID, size, writePositions.data(),
                                   writeIDs.data());
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        // Vector3d force ( Vector3D(counter) + wrap<3,double>(it.vertexCoords()) );
        Vector3d force ( Vector3d::Constant(counter) +
                         Eigen::Map<const Vector3d>(it.vertexCoords()) );
        for (int dim=0; dim<3; dim++) forces[it.vertexID()*3+dim] = force[dim];
        pressures[it.vertexID()] = counter + it.vertexCoords()[0];
      }
      cplInterface.writeBlockVectorData(forcesID, size, writeIDs.data(), forces.data());
      cplInterface.writeBlockScalarData(pressuresID, size, writeIDs.data(), pressures.data());

      cplInterface.getMeshVertices(meshOneID, size, writeIDs.data(),
                                   getWritePositions.data());
      validateWithParams2(math::equals(writePositions, getWritePositions),
                          writePositions, getWritePositions);

      cplInterface.getMeshVertexIDsFromPositions(meshOneID, size, writePositions.data(),
                                                 getWriteIDs.data());
      validateWithParams2(math::equals(writeIDs, getWriteIDs), writeIDs, getWriteIDs);
      //cplInterface.mapWrittenData(meshID);
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          for (int dim=0; dim < 3; dim++){
            int index = it.vertexID()*3+dim;
            readPositions[index] = it.vertexCoords()[dim];
            expectedVelocities[index] = counter + it.vertexCoords()[dim];
          }
          expectedTemperatures[it.vertexID()] = counter + it.vertexCoords()[0];
        }
        cplInterface.resetMesh(meshOneID);
        cplInterface.setMeshVertices(meshOneID, size, readPositions.data(), readIDs.data());
        cplInterface.mapReadDataTo(meshOneID);
        cplInterface.readBlockVectorData(velocitiesID, size, readIDs.data(),
                                         velocities.data());
        cplInterface.readBlockScalarData(temperaturesID, size, readIDs.data(),
                                         temperatures.data());
        validateWithParams2(math::equals(velocities, expectedVelocities),
                            velocities, expectedVelocities);
        validateWithParams2(math::equals(temperatures, expectedTemperatures),
                            temperatures, expectedTemperatures);

        counter += 1.0;
      }
    }
    cplInterface.finalize ();
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    SolverInterface cplInterface ( "SolverTwo", 0, 1 );
    configureSolverInterface ( _pathToTests + "/explicit-mpi-single-non-inc.xml",
                               cplInterface );
    double maxDt = cplInterface.initialize ();
    int squareID = cplInterface.getMeshID("Test-Square");
    int forcesID = cplInterface.getDataID ( "Forces", squareID );
    int pressuresID = cplInterface.getDataID("Pressures", squareID );
    int velocitiesID = cplInterface.getDataID ( "Velocities", squareID );
    int temperaturesID = cplInterface.getDataID("Temperatures", squareID );
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    // SolverTwo does not start the coupled simulation and has, hence,
    // already received the first data to be validated.
    for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
      Vector3d force = Vector3d::Zero();
      double pressure = 0.0;
      cplInterface.readVectorData(forcesID, it.vertexID(), force.data());
      cplInterface.readScalarData(pressuresID, it.vertexID(), pressure);
      validate(math::equals(force, Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords())) );
      validate(math::equals(pressure, counter + it.vertexCoords()[0]));
    }
    counter += 1.0;

    while ( cplInterface.isCouplingOngoing() ) {
      for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
        Vector3d vel ( Vector3d::Constant(counter - 1.0) + Eigen::Map<const Vector3d>(it.vertexCoords()) );
        cplInterface.writeVectorData ( velocitiesID, it.vertexID(), vel.data() );
        double temperature = counter - 1.0 + it.vertexCoords()[0];
        cplInterface.writeScalarData(temperaturesID, it.vertexID(), temperature);
      }
      maxDt = cplInterface.advance ( maxDt );
      if ( cplInterface.isCouplingOngoing() ) {
        for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
          Vector3d force = Vector3d::Zero();
          double pressure = 0.0;
          cplInterface.readVectorData(forcesID, it.vertexID(), force.data());
          cplInterface.readScalarData(pressuresID, it.vertexID(), pressure);
          validate ( math::equals(force, Vector3d::Constant(counter) + Eigen::Map<const Vector3d>(it.vertexCoords())) );
          validate(math::equals(pressure, counter + it.vertexCoords()[0]));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithSolverGeometry ()
{
  TRACE();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );

  mesh::Mesh::resetGeometryIDsGlobally ();

  int timesteps = 0;
  double time = 0;

  if ( utils::Parallel::getProcessRank() == 0 ){
    runSolver ( "SolverOne",
                "/explicit-solvergeometry.xml",
                timesteps, time );
  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    SolverInterface cplInterface ( "SolverTwo", 0, 1 );
    configureSolverInterface (
      _pathToTests + "/explicit-solvergeometry.xml",
      cplInterface );
    validateEquals ( cplInterface.getDimensions(), 3 );
    int meshID = cplInterface.getMeshID ( "SolverGeometry" );
    int i0 = cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0,0.0,0.0).data());
    int i1 = cplInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0,0.0,0.0).data());
    int i2 = cplInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0,1.0,0.0).data());
    int e0 = cplInterface.setMeshEdge(meshID, i0, i1);
    int e1 = cplInterface.setMeshEdge(meshID, i1, i2);
    int e2 = cplInterface.setMeshEdge(meshID, i2, i0);
    cplInterface.setMeshTriangle(meshID, e0, e1, e2);
    double dt = cplInterface.initialize();

    int size = cplInterface.getMeshVertexSize(meshID);
    validateEquals(size, 3);

    while ( cplInterface.isCouplingOngoing() ){
      time += dt;
      dt = cplInterface.advance(dt);
      timesteps++;
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithDisplacingGeometry()
{
  TRACE();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );

  using Eigen::Vector3d;
  using math::equals;
  int timesteps = 0;
  double time = 0;

  Vector3d zero = Vector3d::Zero();
  Vector3d one(1, 0, 0);
  Vector3d two(0, 1, 0);
  Vector3d three(1, 1, 1);
  
  if ( utils::Parallel::getProcessRank() == 0 ) { // SolverOne part
    // SolverOne has a local offset to the geometry
    Vector3d localOffset = Vector3d::Constant(10);
    zero += localOffset;
    one += localOffset;
    two += localOffset;
    three += localOffset;

    SolverInterface cplInterface ( "SolverOne", 0, 1 );
    configureSolverInterface (
      _pathToTests + "/explicit-solvergeometry.xml",
      cplInterface );
    double dt = cplInterface.initialize();

    int meshID = cplInterface.getMeshID("SolverGeometry");
    int size = cplInterface.getMeshVertexSize(meshID);
    validateEquals(size, 4);

    while (cplInterface.isCouplingOngoing()){
      MeshHandle handle = cplInterface.getMeshHandle("SolverGeometry");
      VertexIterator iter = handle.vertices().begin();
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), zero));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), one));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), two));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), three));
      iter++;
      validate(not (iter != handle.vertices().end()) );

      time += dt;
      timesteps++;
      dt = cplInterface.advance(dt);

      // Add displacements (known in this test case) for validation
      zero += Vector3d::Constant(1.0);
      one += Vector3d::Constant(2.0);
      two += Vector3d::Constant(3.0);
      three += Vector3d::Constant(4.0);
    }
    cplInterface.finalize();
  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    SolverInterface cplInterface ( "SolverTwo", 0, 1 );
    configureSolverInterface (
      _pathToTests + "/explicit-solvergeometry.xml",
      cplInterface );
    int meshID = cplInterface.getMeshID("SolverGeometry");
    int i0 = cplInterface.setMeshVertex(meshID, zero.data());
    int i1 = cplInterface.setMeshVertex(meshID, one.data());
    int i2 = cplInterface.setMeshVertex(meshID, two.data());
    int i3 = cplInterface.setMeshVertex(meshID, three.data());
//      cplInterface.setMeshEdge ( meshID, i0, i1 );
//      cplInterface.setMeshEdge ( meshID, i1, i3 );
//      cplInterface.setMeshEdge ( meshID, i3, i2 );
//      cplInterface.setMeshEdge ( meshID, i2, i0 );

    double dt = cplInterface.initialize();
    int displacementsID = cplInterface.getDataID("Displacements", meshID);

    while (cplInterface.isCouplingOngoing()){
      MeshHandle handle = cplInterface.getMeshHandle("SolverGeometry");
      VertexIterator iter = handle.vertices().begin();
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), zero));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), one));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), two));
      iter++;
      validate(equals(Eigen::Map<const Vector3d>(iter.vertexCoords()), three));
      iter++;
      validate(not (iter != handle.vertices().end()));

      // Add displacements
      cplInterface.writeVectorData(displacementsID, i0, Vector3d(1, 1, 1).data());
      cplInterface.writeVectorData(displacementsID, i1, Vector3d(2, 2, 2).data());
      cplInterface.writeVectorData(displacementsID, i2, Vector3d(3, 3, 3).data());
      cplInterface.writeVectorData(displacementsID, i3, Vector3d(4, 4, 4).data());

      // modify coordinates by displacements
      zero += Vector3d::Constant(1);
      one += Vector3d::Constant(2);
      two += Vector3d::Constant(3);
      three += Vector3d::Constant(4);

      time += dt;
      dt = cplInterface.advance(dt);
      timesteps++;
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithDataScaling()
{
  TRACE();
  assertion ( utils::Parallel::getCommunicatorSize() == 2 );
  double dt;
  if ( utils::Parallel::getProcessRank() == 0 ) { // SolverOne part
    SolverInterface cplInterface ( "SolverOne", 0, 1 );
    configureSolverInterface (
        _pathToTests + "/explicit-datascaling.xml",
        cplInterface );
    validateEquals ( cplInterface.getDimensions(), 2 );
    dt = cplInterface.initialize();
    int meshID = cplInterface.getMeshID("Test-Square");
    int velocitiesID = cplInterface.getDataID ( "Velocities", meshID );
    while ( cplInterface.isCouplingOngoing() ){
      MeshHandle handle = cplInterface.getMeshHandle ( "Test-Square" );
      VertexIterator iter = handle.vertices().begin();
      for ( size_t i=0; iter != handle.vertices().end(); i++, iter++ ) {
        Eigen::Vector2d data = Eigen::Vector2d::Constant(i);
        cplInterface.writeVectorData ( velocitiesID, i, data.data() );
      }
      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();

  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    SolverInterface cplInterface ( "SolverTwo", 0, 1 );
    configureSolverInterface (
        _pathToTests + "/explicit-datascaling.xml",
        cplInterface );
    validateEquals ( cplInterface.getDimensions(), 2 );
    dt = cplInterface.initialize();
    int meshID = cplInterface.getMeshID("Test-Square");
    int velocitiesID = cplInterface.getDataID ( "Velocities", meshID );
    while ( cplInterface.isCouplingOngoing() ){
      MeshHandle handle = cplInterface.getMeshHandle ( "Test-Square" );
      VertexIterator iter = handle.vertices().begin();
      for ( size_t i=0; iter != handle.vertices().end(); iter++, i++ ){
        Eigen::Vector2d readData;
        cplInterface.readVectorData ( velocitiesID, i, readData.data() );
        Eigen::Vector2d expectedData = Eigen::Vector2d::Constant(i * 10.0);
        validate ( math::equals(readData, expectedData, 5e-13) );
      }
      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();
  }
}


void SolverInterfaceTest:: testExplicitWithCheckpointingStatMapping()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  using namespace precice::constants;
  int timesteps = 0;
  double time = 0.0;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-writecheckpoint-stat.xml",
                             couplingInterface);
    validateEquals(couplingInterface.getDimensions(), 2);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");

    int meshOneID = couplingInterface.getMeshID("MeshOne");

    double pos[2];
    // Set mesh positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setMeshVertex(meshOneID, pos);

    double dt = couplingInterface.initialize();
    int forcesID = couplingInterface.getDataID("Forces", meshOneID);
    /*int velocitiesID =*/ couplingInterface.getDataID("Velocities", meshOneID);
    validateEquals(solverOne->_meshContexts.size(), 2);

    while (couplingInterface.isCouplingOngoing()){
      validate(solverOne->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverOne->_dataContexts[forcesID];
      validate(dataContext != nullptr);
      mesh::PtrData data = dataContext->fromData;
      validateEquals(forcesID, dataContext->fromData->getID());
      auto& values = data->values();
      values.setConstant(1);
      time += dt;
      dt = couplingInterface.advance(dt);
      couplingInterface.mapReadDataTo(meshOneID);
      timesteps++;
      if (couplingInterface.isActionRequired(actionWriteSimulationCheckpoint())){
        validateEquals(timesteps, 10);
        couplingInterface.fulfilledAction(actionWriteSimulationCheckpoint());
      }
    }
    couplingInterface.finalize();
    validateEquals(timesteps, 10);
    validateEquals(time, 10.0);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-writecheckpoint-stat.xml",
                             couplingInterface);
    validateEquals(couplingInterface.getDimensions(), 2);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    double dt = couplingInterface.initialize();
    int squareID = couplingInterface.getMeshID("Test-Square");
    int dataID = couplingInterface.getDataID("Velocities", squareID);
    validateEquals(solverTwo->_meshContexts.size(), 2);
    while (couplingInterface.isCouplingOngoing()){
      validate(solverTwo->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverTwo->_dataContexts[dataID];
      validate(dataContext != nullptr);
      mesh::PtrData fromData = dataContext->fromData;
      mesh::PtrData toData = dataContext->toData;
      // no mapping, so fromData and toData should be the same data
      validateEquals(fromData->getID(), toData->getID());
      validateEquals(dataID, fromData->getID());
      auto& values = fromData->values();
      values.setConstant(2);
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
      if (couplingInterface.isActionRequired(actionWriteSimulationCheckpoint())){
        validateEquals(timesteps, 10);
        couplingInterface.fulfilledAction(actionWriteSimulationCheckpoint());
      }
    }
    couplingInterface.finalize();
    validateEquals(timesteps, 10);
    validateEquals(time, 10.0);
  }

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-readcheckpoint-stat.xml",
                             couplingInterface);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");

    int meshOneID = couplingInterface.getMeshID("MeshOne");

    double pos[2];
    // Set mesh positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setMeshVertex(meshOneID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setMeshVertex(meshOneID, pos);

    double dt = couplingInterface.initialize();
    /*int forcesID = */ couplingInterface.getDataID("Forces", meshOneID);
    validateEquals(solverOne->_meshContexts.size(), 2);
    mesh::PtrMesh mesh = solverOne->_meshContexts[0]->mesh;
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());

    while (couplingInterface.isCouplingOngoing()){
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();

    validateEquals(timesteps, 20);
    validateEquals(time, 20.0);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-readcheckpoint-stat.xml",
                             couplingInterface);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    double dt = couplingInterface.initialize();
    int squareID = couplingInterface.getMeshID("Test-Square");
    /*int dataID = */ couplingInterface.getDataID("Velocities", squareID);
    validateEquals(solverTwo->_meshContexts.size(), 2);
    mesh::PtrMesh mesh = solverTwo->_meshContexts[0]->mesh;
    Eigen::Vector2d integral(0, 0);
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());
    while (couplingInterface.isCouplingOngoing()){
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();
    validateEquals(timesteps, 20);
    validateEquals(time, 20.0);
  }
}

void SolverInterfaceTest:: testImplicit()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  double state = 0.0;
  double checkpoint = 0.0;
  int iterationCount = 0;
  double initialStateChange = 5.0;
  double stateChange = initialStateChange;
  int computedTimesteps = 0;
  using namespace precice::constants;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit.xml", couplingInterface);
    double maxDt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
      }
      iterationCount++;
      stateChange = initialStateChange / (double)iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimestepComplete()){
        computedTimesteps ++;
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 4);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit.xml", couplingInterface);
    double maxDt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
        iterationCount++;
      }
      stateChange = initialStateChange / (double)iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimestepComplete()){
        computedTimesteps++;
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 4);
  }
}


void SolverInterfaceTest:: testImplicitWithCheckpointingMappingStat()
{
  TRACE()
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  double state = 0.0;
  double checkpoint = 0.0;
  int iterationCount = 0;
  double initialStateChange = 5.0;
  double stateChange = initialStateChange;
  int computedTimesteps = 0;
  using namespace precice::constants;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit-writecheckpoint-stat.xml", couplingInterface);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    int meshID = couplingInterface.getMeshID("Square");
    int forcesID = couplingInterface.getDataID("Forces", meshID);
    int velocitiesID = couplingInterface.getDataID("Velocities", meshID);

    double maxDt = couplingInterface.initialize();

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
      }
      validate(solverOne->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverOne->_dataContexts[forcesID];
      validate(dataContext != nullptr);
      mesh::PtrData localForces = dataContext->fromData;
      validateEquals(forcesID, localForces->getID());
      auto& forceValues = localForces->values();
      forceValues.setConstant(1);
      iterationCount++;
      maxDt = couplingInterface.advance(maxDt);

      dataContext = solverOne->_dataContexts[velocitiesID];
      validate(dataContext != nullptr);
      mesh::PtrData localVelocities = dataContext->toData;
      validateEquals(velocitiesID, localVelocities->getID());
      auto& velocityValues = localVelocities->values();
      Eigen::Vector2d integral;
      math::sumSubvectors(velocityValues, integral);
      validate(math::equals(integral, Eigen::Vector2d(8.0, 8.0)));

      if (couplingInterface.isTimestepComplete()){
        computedTimesteps ++;
      }
      if (couplingInterface.isActionRequired(actionWriteSimulationCheckpoint())){
        validateEquals(computedTimesteps, 4);
        couplingInterface.fulfilledAction(actionWriteSimulationCheckpoint());
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 4);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit-writecheckpoint-stat.xml", couplingInterface);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    int meshID = couplingInterface.getMeshID("Square");
    int forcesID = couplingInterface.getDataID("Forces", meshID);
    int velocitiesID = couplingInterface.getDataID("Velocities", meshID);
    double maxDt = couplingInterface.initialize();

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
        iterationCount++;
      }
      validate(solverTwo->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverTwo->_dataContexts[velocitiesID];
      validate(dataContext != nullptr);
      mesh::PtrData velocities = dataContext->fromData;
      validateEquals(velocitiesID, velocities->getID());
      auto& values = velocities->values();
      values = Eigen::VectorXd::Constant(values.size(), 2.0);
      
      maxDt = couplingInterface.advance(maxDt);

      dataContext = solverTwo->_dataContexts[forcesID];
      validate(dataContext != nullptr);
      mesh::PtrData forces = dataContext->toData;
      //validateEquals(velocitiesID, localVelocities->getID());
      auto& forceValues = forces->values();
      Eigen::Vector2d integral;
      math::sumSubvectors(forceValues, integral);
      validate(math::equals(integral, Eigen::Vector2d(4.0, 4.0)));

      if (couplingInterface.isTimestepComplete()){
        computedTimesteps++;
      }
      if (couplingInterface.isActionRequired(actionWriteSimulationCheckpoint())){
        validateEquals(computedTimesteps, 4);
        couplingInterface.fulfilledAction(actionWriteSimulationCheckpoint());
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 4);
  }

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit-readcheckpoint-stat.xml", couplingInterface);
    int meshID = couplingInterface.getMeshID("Square");
    /*int forcesID =*/ couplingInterface.getDataID("Forces", meshID);
    int velocitiesID = couplingInterface.getDataID("Velocities", meshID);
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());

    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    impl::DataContext* dataContext = solverOne->_dataContexts[velocitiesID];
    validate(dataContext != nullptr);
    mesh::PtrData localVelocities = dataContext->toData;
    validateEquals(velocitiesID, localVelocities->getID());
    auto& velocityValues = localVelocities->values();
    Eigen::Vector2d integral;
    math::sumSubvectors(velocityValues, integral);
    validate(math::equals(integral, Eigen::Vector2d(8.0, 8.0)));

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
      }
      iterationCount ++;
      stateChange = initialStateChange / (double)iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimestepComplete()){
        computedTimesteps ++;
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 6);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/implicit-readcheckpoint-stat.xml", couplingInterface);
    int meshID = couplingInterface.getMeshID("Square");
    int forcesID = couplingInterface.getDataID("Forces", meshID);
    /*int velocitiesID =*/ couplingInterface.getDataID("Velocities", meshID);
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());

    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    impl::DataContext* dataContext = solverTwo->_dataContexts[forcesID];
    validate(dataContext != nullptr);
    mesh::PtrData forces = dataContext->toData;
    validateEquals(forcesID, forces->getID());
    auto& forceValues = forces->values();
    Eigen::Vector2d integral;
    math::sumSubvectors(forceValues, integral);
    validate(math::equals(integral, Eigen::Vector2d(4.0, 4.0)));

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionReadIterationCheckpoint());
        state = checkpoint;
        iterationCount ++;
      }
      stateChange = initialStateChange / (double)iterationCount;
      state += stateChange;
      maxDt = couplingInterface.advance(maxDt);
      if (couplingInterface.isTimestepComplete()){
        computedTimesteps++;
      }
    }
    couplingInterface.finalize();
    validateEquals(computedTimesteps, 6);
  }
}

void SolverInterfaceTest:: runSolver
(
   const std::string& solverName,
   const std::string& configurationFileName,
   int&               timestepsComputed,
   double&            timeComputed )
{
  TRACE(solverName, configurationFileName);
  timestepsComputed = 0;
  timeComputed = 0.0;
  SolverInterface couplingInterface(solverName, 0, 1);
  configureSolverInterface(_pathToTests + configurationFileName, couplingInterface);
  validateEquals(couplingInterface.getDimensions(), 3);
  double dt = couplingInterface.initialize();
  while (couplingInterface.isCouplingOngoing()){
    timeComputed += dt;
    dt = couplingInterface.advance(dt);
    timestepsComputed++;
  }
  couplingInterface.finalize();
}

void SolverInterfaceTest:: testStationaryMappingWithSolverMesh()
{
  TRACE();
  std::string config2D = _pathToTests + "mapping-without-geo-2D.xml";
  std::string config3D = _pathToTests + "mapping-without-geo-3D.xml";
  int rank = utils::Parallel::getProcessRank();
  assertion((rank == 0) || (rank == 1), rank);
  std::string solverName = rank == 0 ? "SolverA" : "SolverB";
  std::string meshForcesA = "MeshForcesA";
  std::string meshDisplA = "MeshDisplacementsA";
  std::string meshForcesB = "MeshForcesB";
  std::string meshDisplB = "MeshDisplacementsB";
  std::string dataForces = constants::dataForces();
  std::string dataDispl = constants::dataDisplacements();
  using math::equals;

  for (int dim=2; dim < 3; dim++){
    DEBUG("Running " << dim << "D test");
    SolverInterface interface(solverName, 0, 1);
    if (dim == 2){
      configureSolverInterface(config2D, interface);
    }
    else {
      configureSolverInterface(config3D, interface);
    }
    validateEquals(interface.getDimensions(), dim);

    std::vector<Eigen::VectorXd> positions;
    Eigen::VectorXd position(dim);
    if (dim == 2){
      position << 0.0, 0.0;
      positions.push_back(position);
      position << 1.0, 0.0;
      positions.push_back(position);
      position << 1.0, 1.0;
      positions.push_back(position);
      position << 0.0, 1.0;
      positions.push_back(position);
    }
    else {
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

    if (rank == 0){
      int meshForcesID = interface.getMeshID(meshForcesA);
      int meshDisplID = interface.getMeshID(meshDisplA);
      int dataForcesID = interface.getDataID(dataForces, meshForcesID);
      int dataDisplID = interface.getDataID(dataDispl, meshDisplID);

      // Set solver mesh positions for reading and writing data with mappings
      for (size_t i=0; i < size; i++){
        position = positions[i].array() + 0.1;
        interface.setMeshVertex(meshForcesID, position.data());
        position = positions[i].array() + 0.6;
        interface.setMeshVertex(meshDisplID, position.data());
      }
      double maxDt = interface.initialize();

      validate(interface.isWriteDataRequired(maxDt));
      validate(not interface.isReadDataAvailable());
      Eigen::VectorXd force = Eigen::VectorXd::Constant(dim, 1);
      Eigen::VectorXd displ = Eigen::VectorXd::Constant(dim, 0);
      for (size_t i=0; i < size; i++){
        interface.writeVectorData(dataForcesID, i, force.data());
      }
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      interface.mapReadDataTo(meshDisplID);
      //INFO("1: mapped data: " << interface._impl->_accessor->dataContext(dataDisplID).data->values());
      force.array() += 1.0;
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, displ.data());
        validateNumericalEquals(displ[0], positions[i][0] + 0.1);
        interface.writeVectorData(dataForcesID, i, force.data());
      }
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      interface.mapReadDataTo(meshDisplID);
      //INFO("2: mapped data: " << interface._impl->_accessor->dataContext(dataDisplID).data->values());
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, displ.data());
        validateNumericalEquals(displ[0], 2.0*(positions[i][0] + 0.1));
      }
      interface.finalize();
    }
    else {
      assertion(rank == 1, rank);
      int meshForcesID = interface.getMeshID(meshForcesB);
      int meshDisplID = interface.getMeshID(meshDisplB);
      int dataForcesID = interface.getDataID(dataForces, meshForcesID);
      int dataDisplID = interface.getDataID(dataDispl, meshDisplID);

      // Set solver mesh positions provided to SolverA for data mapping
      for (size_t i=0; i < size; i++){
        interface.setMeshVertex(meshForcesID, positions[i].data());
        position = positions[i].array() + 0.5;
        interface.setMeshVertex(meshDisplID, position.data());
      }
      double maxDt = interface.initialize();

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      Eigen::VectorXd force = Eigen::VectorXd::Zero(dim);
      Eigen::VectorXd totalForce = Eigen::VectorXd::Zero(dim);
      Eigen::VectorXd displ = Eigen::VectorXd::Zero(dim);
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataForcesID, i, force.data());
        totalForce += force;
        displ.setConstant(positions[i][0]);
        interface.writeVectorData(dataDisplID, i, displ.data());
      }
      Eigen::VectorXd expected = Eigen::VectorXd::Constant(dim, size);
      validateWithParams2(math::equals(totalForce,expected), totalForce, expected);
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      totalForce.setConstant(0);
      for (size_t i=0; i < positions.size(); i++){
        interface.readVectorData(dataForcesID, i, force.data());
        totalForce += force;
        displ.setConstant(2.0 * positions[i][0]);
        interface.writeVectorData(dataDisplID, i, displ.data());
      }
      expected.setConstant(2.0 * (double)size);
      validateWithParams1(math::equals(totalForce,expected), totalForce);
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(not interface.isReadDataAvailable()); //second participant has no new data after last advance
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, force.data());
      }
      interface.finalize();
    }
  }
}

void SolverInterfaceTest:: testDistributedCommunications()
{
  TRACE();

  assertion(utils::Parallel::getCommunicatorSize() == 4);

  std::vector<std::string> fileNames({
      "point-to-point-sockets.xml",
      "point-to-point-mpi.xml",
      "gather-scatter-mpi.xml"});

  for (auto fileName : fileNames) {
    mesh::Mesh::resetGeometryIDsGlobally();

    std::string solverName;
    int rank = -1, size = -1;
    std::string meshName;
    int i1 = -1 ,i2 = -1; //indices for data and positions

    std::vector<Eigen::VectorXd> positions;
    std::vector<Eigen::VectorXd> data;
    std::vector<Eigen::VectorXd> expectedData;

    Eigen::Vector3d position;
    Eigen::Vector3d datum;

    for( int i=0; i<4; i++){
      position[0] = i*1.0;
      position[1] = 0.0;
      position[2] = 0.0;
      positions.push_back(position);
      datum[0] = i*1.0;
      datum[1] = i*1.0;
      datum[2] = 0.0;
      data.push_back(datum);
      datum[0] = i*2.0+1.0;
      datum[1] = i*2.0+1.0;
      datum[2] = 1.0;
      expectedData.push_back(datum);
    }

    if (utils::Parallel::getProcessRank() == 0){
      solverName = "Fluid";
      rank = 0;
      size = 2;
      meshName = "FluidMesh";
      i1 = 0;
      i2 = 2;
    }
    else if(utils::Parallel::getProcessRank() == 1){
      solverName = "Fluid";
      rank = 1;
      size = 2;
      meshName = "FluidMesh";
      i1 = 2;
      i2 = 4;
    }
    else if(utils::Parallel::getProcessRank() == 2){
      solverName = "Structure";
      rank = 0;
      size = 2;
      meshName = "StructureMesh";
      i1 = 0;
      i2 = 1;
    }
    else if(utils::Parallel::getProcessRank() == 3){
      solverName = "Structure";
      rank = 1;
      size = 2;
      meshName = "StructureMesh";
      i1 = 1;
      i2 = 4;
    }

    SolverInterface precice(solverName, rank, size);
    configureSolverInterface(_pathToTests + fileName, precice);
    int meshID = precice.getMeshID(meshName);
    int forcesID = precice.getDataID("Forces", meshID);
    int velocID = precice.getDataID("Velocities", meshID);

    std::vector<int> vertexIDs;
    for(int i=i1; i<i2; i++){
      int vertexID = precice.setMeshVertex(meshID, positions[i].data());
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    if (utils::Parallel::getProcessRank() <= 1){ //Fluid
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.writeVectorData(forcesID, vertexIDs[i], data[i+i1].data());
      }
    }
    else if (utils::Parallel::getProcessRank() >= 2){ //Structure
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.readVectorData(forcesID, vertexIDs[i], data[i].data());
        data[i] = (data[i]*2).array() + 1.0;
        precice.writeVectorData(velocID, vertexIDs[i], data[i].data());
      }
    }

    precice.advance(1.0);

    if (utils::Parallel::getProcessRank() <= 1){ //Fluid
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.readVectorData(velocID, vertexIDs[i], data[i+i1].data());
        for (size_t d=0; d<3; d++){
          validateNumericalEquals(expectedData[i+i1][d],data[i+i1][d]);
        }
      }
    }

    precice.finalize();
  }
}

void SolverInterfaceTest:: testBug()
{
  TRACE();
  using Eigen::Vector3d;
  std::string config = _pathToTests + "bug.xml";

  int slices = 5;
  std::vector<Vector3d> coords;
  for (int i=0; i < slices; i++){
    double z = (double)i * 1.0;
    coords.push_back( Vector3d( 1.0,  0.0, z) );
    coords.push_back( Vector3d( 0.0,  1.0, z) );
    coords.push_back( Vector3d(-1.0,  0.0, z) );
    coords.push_back( Vector3d( 0.0, -1.0, z) );
  }
      
  int rank = utils::Parallel::getProcessRank();
  assertion((rank == 0) || (rank == 1), rank);
  std::string solverName = rank == 0 ? "Flite" : "Calculix";
  if (solverName == std::string("Flite")){
    SolverInterface precice("Flite", 0, 1);
    configureSolverInterface(config, precice);
    int meshID = precice.getMeshID("FliteNodes");
    int forcesID = precice.getDataID(precice::constants::dataForces(), meshID);
    int displacementsID = precice.getDataID(precice::constants::dataDisplacements(), meshID);
    int oldDisplacementsID = precice.getDataID("OldDisplacements", meshID);
    validateEquals(precice.getDimensions(), 3);
    for (Vector3d& coord : coords){
      precice.setMeshVertex(meshID, coord.data());
    }
    double maxDt = precice.initialize();
    double dt = 1.0e-5 / 15.0; // Flite took 15 subcycling steps
    while (precice.isCouplingOngoing()){
      dt = dt < maxDt ? dt : maxDt;
      for(int i=0; i < (int)coords.size(); i++){
        double force[3] = { 1.0, 2.0, 3.0 };
        precice.writeVectorData(forcesID, i, force);
      }
      maxDt = precice.advance(dt);
      precice.mapReadDataTo(meshID);
      for(int i=0; i < (int)coords.size(); i++){
        double displacement[3];
        double oldDisplacement[3];
        precice.readVectorData(displacementsID, i, displacement);
        precice.readVectorData(oldDisplacementsID, i, oldDisplacement);
      }
    }
    precice.finalize();
  }
  else {
    assertion(solverName == std::string("Calculix"), solverName);
    SolverInterface precice("Calculix", 0, 1);
    configureSolverInterface(config, precice);
    int meshID = precice.getMeshID("CalculixNodes");
    for (Vector3d& coord : coords){
      precice.setMeshVertex(meshID, coord.data());
    }
    for(int i=0; i < slices-1; i++){
      // Build cylinder/channel geometry
      precice.setMeshTriangleWithEdges(meshID, i*4, (i*4)+1, (i+1)*4);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4, (i*4)+1, ((i+1)*4)+1);
      precice.setMeshTriangleWithEdges(meshID, i*4+1, (i*4)+2, (i+1)*4+1);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+1, (i*4)+2, ((i+1)*4)+2);
      precice.setMeshTriangleWithEdges(meshID, i*4+2, (i*4)+3, (i+1)*4+2);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+2, (i*4)+3, ((i+1)*4)+3);
      precice.setMeshTriangleWithEdges(meshID, i*4+3, (i*4), (i+1)*4+3);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+3, i*4, (i+1)*4);
    }
    double dt = precice.initialize();
    while(precice.isCouplingOngoing()){
      precice.advance(dt);
    }
    precice.finalize();
  }
}

void SolverInterfaceTest:: testThreeSolvers()
{
  TRACE();
  std::string configFilename(_pathToTests + "three-solver-explicit-explicit.xml");
  std::vector<int> expectedCallsOfAdvance = {10, 10, 10};
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-implicit-implicit.xml";
  expectedCallsOfAdvance = {30, 30, 20};
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-implicit-explicit.xml";
  expectedCallsOfAdvance = {30, 30, 10};
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-explicit-implicit.xml";
  expectedCallsOfAdvance = {30, 10, 30};
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-parallel.xml";
  expectedCallsOfAdvance = {30, 30, 10};
  runThreeSolvers(configFilename, expectedCallsOfAdvance);
}

void SolverInterfaceTest:: runThreeSolvers
(
  const std::string&      configFilename,
  const std::vector<int>& expectedCallsOfAdvance )
{
  TRACE(configFilename, expectedCallsOfAdvance);

  int rank = utils::Parallel::getProcessRank();
  assertion((rank == 0) || (rank == 1) || (rank == 2), rank);

  std::string writeIterCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterCheckpoint(constants::actionReadIterationCheckpoint());
  std::string writeInitData(constants::actionWriteInitialData());

  std::string solverName;
  if (rank == 0) solverName = std::string("SolverOne");
  else if (rank == 1) solverName = std::string("SolverTwo");
  else solverName = std::string("SolverThree");
  int callsOfAdvance = 0;

  if (solverName == std::string("SolverOne")){
    SolverInterface precice(solverName, 0, 1);
    configureSolverInterface(configFilename, precice);
    int meshID = precice.getMeshID("Mesh");
    //int dataID = precice.getDataID("Data");
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)){
      precice.fulfilledAction(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()){
      //precice.writeVectorData(dataID, 0, raw(Vector2D(1.0, 2.0)));
      if (precice.isActionRequired(writeIterCheckpoint)){
        precice.fulfilledAction(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)){
        precice.fulfilledAction(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    validateEquals(callsOfAdvance, expectedCallsOfAdvance[0]);
  }
  else if (solverName == std::string("SolverTwo")){
    SolverInterface precice(solverName, 0, 1);
    configureSolverInterface(configFilename, precice);
    //int dataID = precice.getDataID("Data");
    //precice.setReadPosition(meshID, raw(utils::Vector2D(0.0, 0.0))); //no use here
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)){
      precice.fulfilledAction(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()){
      //Vector2D data;
      //precice.readVectorData(dataID, 0, raw(data));
      if (precice.isActionRequired(writeIterCheckpoint)){
        precice.fulfilledAction(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)){
        precice.fulfilledAction(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    validateEquals(callsOfAdvance, expectedCallsOfAdvance[1]);
  }
  else {
    assertion(solverName == std::string("SolverThree"), solverName);
    SolverInterface precice(solverName, 0, 1);
    configureSolverInterface(configFilename, precice);
    //int dataID = precice.getDataID("Data");
    //precice.setReadPosition(meshID, raw(utils::Vector2D(0.0, 0.0))); //no use here
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)){
      precice.fulfilledAction(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()){
      //Vector2D data;
      //precice.readVectorData(dataID, 0, raw(data));
      if (precice.isActionRequired(writeIterCheckpoint)){
        precice.fulfilledAction(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)){
        precice.fulfilledAction(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    validateEquals(callsOfAdvance, expectedCallsOfAdvance[2]);
  }
}

void SolverInterfaceTest:: testMultiCoupling()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  mesh::Mesh::resetGeometryIDsGlobally();

  std::vector<Eigen::Vector2d> positions;
  Eigen::Vector2d position;
  position << 0.0, 0.0;
  positions.push_back(position);
  position << 1.0, 0.0;
  positions.push_back(position);
  position << 1.0, 1.0;
  positions.push_back(position);
  position << 0.0, 1.0;
  positions.push_back(position);

  std::vector<Eigen::Vector2d> datas;
  Eigen::Vector2d data;
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

  if (utils::Parallel::getProcessRank() < 3){
    int meshID = -1;
    int dataWriteID = -1;
    int dataReadID = -1;

    std::string participant = "";

    if (utils::Parallel::getProcessRank() == 0){
      participant = "SOLIDZ1";
    }
    else if (utils::Parallel::getProcessRank() == 1){
      participant = "SOLIDZ2";
    }
    else if (utils::Parallel::getProcessRank() == 2){
      participant = "SOLIDZ3";
    }

    SolverInterface precice(participant, 0, 1);
    configureSolverInterface(_pathToTests + "/multi.xml", precice);
    validateEquals(precice.getDimensions(),2);

    if (utils::Parallel::getProcessRank() == 0){
      meshID = precice.getMeshID("SOLIDZ_Mesh1");
      dataWriteID = precice.getDataID("Displacements1", meshID);
      dataReadID = precice.getDataID("Forces1", meshID);
    }
    else if (utils::Parallel::getProcessRank() == 1){
      meshID = precice.getMeshID("SOLIDZ_Mesh2");
      dataWriteID = precice.getDataID("Displacements2", meshID);
      dataReadID = precice.getDataID("Forces2", meshID);
    }
    else if (utils::Parallel::getProcessRank() == 2){
      meshID = precice.getMeshID("SOLIDZ_Mesh3");
      dataWriteID = precice.getDataID("Displacements3", meshID);
      dataReadID = precice.getDataID("Forces3", meshID);
    }

    std::vector<int> vertexIDs;
    int vertexID = -1;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID, positions[i].data());
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i=0; i < 4; i++){
      precice.writeVectorData(dataWriteID, vertexIDs[i], datas[i].data());
    }

    if (precice.isActionRequired(writeIterCheckpoint)){
      precice.fulfilledAction(writeIterCheckpoint);
    }
    precice.advance(0.0001);
    if (precice.isActionRequired(readIterCheckpoint)){
      precice.fulfilledAction(readIterCheckpoint);
    }

    for (size_t i=0; i < 4; i++){
      precice.readVectorData(dataReadID, vertexIDs[i], datas[i].data());
    }

    validateNumericalEquals(datas[0][0],1.00000000000000002082e-03);
    validateNumericalEquals(datas[0][1],1.00000000000000002082e-03);
    validateNumericalEquals(datas[1][0],0.00000000000000000000e+00);
    validateNumericalEquals(datas[1][1],1.00000000000000002082e-03);
    validateNumericalEquals(datas[2][0],3.00000000000000006245e-03);
    validateNumericalEquals(datas[2][1],3.00000000000000006245e-03);
    validateNumericalEquals(datas[3][0],4.00000000000000008327e-03);
    validateNumericalEquals(datas[3][1],5.00000000000000010408e-03);

    precice.finalize();

  }
  else {
    assertion(utils::Parallel::getProcessRank() == 3);
    SolverInterface precice("NASTIN", 0, 1);
    configureSolverInterface(_pathToTests + "/multi.xml", precice);
    validateEquals(precice.getDimensions(),2);
    int meshID1 = precice.getMeshID("NASTIN_Mesh1");
    int meshID2 = precice.getMeshID("NASTIN_Mesh2");
    int meshID3 = precice.getMeshID("NASTIN_Mesh3");
    int dataWriteID1 = precice.getDataID("Forces1", meshID1);
    int dataWriteID2 = precice.getDataID("Forces2", meshID2);
    int dataWriteID3 = precice.getDataID("Forces3", meshID3);

    std::vector<int> vertexIDs1;
    int vertexID = -1;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID1, positions[i].data());
      vertexIDs1.push_back(vertexID);
    }
    std::vector<int> vertexIDs2;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID2, positions[i].data());
      vertexIDs2.push_back(vertexID);
    }
    std::vector<int> vertexIDs3;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID3, positions[i].data());
      vertexIDs3.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i=0; i < 4; i++){
      precice.writeVectorData(dataWriteID1, vertexIDs1[i], datas[i].data());
      precice.writeVectorData(dataWriteID2, vertexIDs2[i], datas[i].data());
      precice.writeVectorData(dataWriteID3, vertexIDs3[i], datas[i].data());
    }

    if (precice.isActionRequired(writeIterCheckpoint)){
      precice.fulfilledAction(writeIterCheckpoint);
    }
    precice.advance(0.0001);
    if (precice.isActionRequired(readIterCheckpoint)){
      precice.fulfilledAction(readIterCheckpoint);
    }

    precice.finalize();

  }

}

void SolverInterfaceTest:: testNASTINMeshRestart()
{
  TRACE();
  assertion(utils::Parallel::getCommunicatorSize() == 2);

  std::vector<std::string> restartFiles;
  restartFiles.push_back("precice_checkpoint_NASTIN_NASTIN_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_NASTIN_SOLIDZ_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_NASTIN_simstate.txt");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_cplscheme.txt");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_SOLIDZ_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_simstate.txt");

  for (std::string& restartFile : restartFiles){
    std::ifstream  src(_pathToTests + restartFile, std::ifstream::in);
    std::ofstream  dst(restartFile, std::ifstream::out);
    dst << src.rdbuf();
  }

  std::string readSimCheckpoint(constants::actionReadSimulationCheckpoint());
  std::string writeItCheckpoint(constants::actionWriteIterationCheckpoint());

  mesh::Mesh::resetGeometryIDsGlobally();
  int meshSize = 27;

  std::vector<double> positions(meshSize*2);

  int meshID = -1;

  std::string participant = "";

  if (utils::Parallel::getProcessRank() == 0){
    participant = "NASTIN";
  }
  else if (utils::Parallel::getProcessRank() == 1){
    participant = "SOLIDZ";
  }

  SolverInterface precice(participant, 0, 1);
  configureSolverInterface(_pathToTests + "/nastin-restart-config.xml", precice);
  validateEquals(precice.getDimensions(),2);

  if (utils::Parallel::getProcessRank() == 0){
    meshID = precice.getMeshID("NASTIN_Mesh");
    positions[0] = 3.0000000000000000;
    positions[1] = 0.59999999999999998;
    positions[2] = 4.0000000000000000;
    positions[3] = 1.3999999999999999;
    positions[4] = 4.0000000000000000;
    positions[5] = 1.6000000000000001;
    positions[6] = 3.0000000000000000;
    positions[7] = 1.8000000000000000;
    positions[8] = 3.0000000000000000;
    positions[9] = 1.6000000000000001;
    positions[10] = 3.0000000000000000;
    positions[11] = 2.0000000000000000;
    positions[12] = 4.0000000000000000,
    positions[13] = 0.80000000000000004;
    positions[14] = 4.0000000000000000;
    positions[15] = 1.0000000000000000;
    positions[16] = 3.0000000000000000;
    positions[17] = 0.0000000000000000;
    positions[18] = 3.0000000000000000;
    positions[19] = 0.20000000000000001;
    positions[20] = 4.0000000000000000;
    positions[21] = 0.40000000000000002;
    positions[22] = 4.0000000000000000;
    positions[23] = 0.59999999999999998;
    positions[24] = 4.0000000000000000;
    positions[25] = 1.2000000000000000;
    positions[26] = 3.0000000000000000;
    positions[27] = 0.40000000000000002;
    positions[28] = 4.0000000000000000;
    positions[29] = 1.8000000000000000;
    positions[30] = 4.0000000000000000;
    positions[31] = 0.20000000000000001;
    positions[32] = 3.0000000000000000;
    positions[33] = 1.3999999999999999;
    positions[34] = 3.0000000000000000;
    positions[35] = 1.2000000000000000;
    positions[36] = 3.0000000000000000;
    positions[37] = 1.0000000000000000;
    positions[38] = 3.0000000000000000;
    positions[39] = 0.80000000000000004;
    positions[40] = 3.5000000000000000;
    positions[41] = 2.0000000000000000;
    positions[42] = 4.0000000000000000;
    positions[43] = 2.0000000000000000;
    positions[44] = 3.8332999999999999;
    positions[45] = 2.0000000000000000;
    positions[46] = 3.6667000000000001;
    positions[47] = 2.0000000000000000;
    positions[48] = 3.1667000000000001;
    positions[49] = 2.0000000000000000;
    positions[50] = 3.3332999999999999;
    positions[51] = 2.0000000000000000;
    positions[52] = 4.0000000000000000;
    positions[53] = 0.0000000000000000;
  }
  else if (utils::Parallel::getProcessRank() == 1){
    meshID = precice.getMeshID("SOLIDZ_Mesh");
    positions[0] = 4.0000000000000000;
    positions[1] = 0.40000000000000002;
    positions[2] = 4.0000000000000000;
    positions[3] = 0.59999999999999998;
    positions[4] = 4.0000000000000000;
    positions[5] = 1.3999999999999999;
    positions[6] = 4.0000000000000000;
    positions[7] = 1.6000000000000001;
    positions[8] = 3.0000000000000000;
    positions[9] = 1.6000000000000001;
    positions[10] = 3.0000000000000000;
    positions[11] = 1.3999999999999999;
    positions[12] = 3.0000000000000000;
    positions[13] = 0.59999999999999998;
    positions[14] = 3.0000000000000000;
    positions[15] = 0.40000000000000002;
    positions[16] = 3.0000000000000000;
    positions[17] = 1.0000000000000000;
    positions[18] = 4.0000000000000000;
    positions[19] = 1.0000000000000000;
    positions[20] = 3.0000000000000000;
    positions[21] = 0.80000000000000004;
    positions[22] = 4.0000000000000000;
    positions[23] = 1.2000000000000000;
    positions[24] = 4.0000000000000000;
    positions[25] = 0.80000000000000004;
    positions[26] = 3.3332999999999999;
    positions[27] = 2.0000000000000000;
    positions[28] = 3.1667000000000001;
    positions[29] = 2.0000000000000000;
    positions[30] = 3.8332999999999999;
    positions[31] = 2.0000000000000000;
    positions[32] = 3.6667000000000001;
    positions[33] = 2.0000000000000000;
    positions[34] = 3.0000000000000000;
    positions[35] = 1.2000000000000000;
    positions[36] = 3.0000000000000000;
    positions[37] = 0.0000000000000000;
    positions[38] = 3.0000000000000000;
    positions[39] = 0.20000000000000001;
    positions[40] = 4.0000000000000000;
    positions[41] = 1.8000000000000000;
    positions[42] = 4.0000000000000000;
    positions[43] = 2.0000000000000000;
    positions[44] = 3.0000000000000000;
    positions[45] = 1.8000000000000000;
    positions[46] = 3.0000000000000000;
    positions[47] = 2.0000000000000000;
    positions[48] = 4.0000000000000000;
    positions[49] = 0.0000000000000000;
    positions[50] = 4.0000000000000000;
    positions[51] = 0.20000000000000001;
    positions[52] = 3.5000000000000000;
    positions[53] = 2.0000000000000000;
  }


  if (precice.isActionRequired(readSimCheckpoint)){
    precice.fulfilledAction(readSimCheckpoint);
  }
  int vertexIDs[meshSize];
  precice.setMeshVertices(meshID, meshSize, positions.data(), vertexIDs);
  precice.initialize();

  if (precice.isActionRequired(writeItCheckpoint)){
    precice.fulfilledAction(writeItCheckpoint);
  }
  precice.finalize();
}

void SolverInterfaceTest:: testPinelliCoupled()
{
  TRACE();

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface interface("EOF", 0, 1);

    {
      std::string filename = "computeForce.py";
      std::ifstream  src((_pathToTests + filename).c_str(), std::ifstream::in);
      std::ofstream  dst(filename.c_str(), std::ifstream::out);
      dst << src.rdbuf();
    }
    configureSolverInterface(_pathToTests + "testPinelliCoupled.xml", interface);
    validateEquals(interface.getDimensions(), 2);
    int meshIdEOF = interface.getMeshID("EOFMesh");
    int dataIdEOFVeloc = interface.getDataID("Velocities",meshIdEOF);
    int dataIdEOFForces = interface.getDataID("Forces",meshIdEOF);

    using Eigen::Vector2d;
    std::vector<int> vertexIdsEOF;

    Vector2d position ( 0.15, 0.15 );
    int vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.15, 0.2;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.15, 0.25;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.2, 0.15;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.2, 0.25;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.25, 0.15;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.25, 0.2;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    position << 0.25, 0.25;
    vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
    vertexIdsEOF.push_back(vertexId);

    interface.initialize();

    for (int vertexID : vertexIdsEOF){
      double data[2] = {2.0,2.0};
      interface.writeVectorData(dataIdEOFVeloc,vertexID,data);
    }

    interface.advance(0.1);

    Vector2d totalForce ( 0.0, 0.0 );

    for (int vertexID : vertexIdsEOF){
      double data[2] = {0.0,0.0};
      interface.readVectorData(dataIdEOFForces,vertexID,data);
      totalForce[0] += data[0];
      totalForce[1] += data[1];
    }

    validateNumericalEquals(totalForce[0]+totalForce[1], 11.0);

    interface.finalize();
  }
  else{
    assertion(utils::Parallel::getProcessRank() == 1);
    SolverInterface interface("Structure", 0, 1);
    configureSolverInterface(_pathToTests + "testPinelliCoupled.xml", interface);
    validateEquals(interface.getDimensions(), 2);
    int meshID = interface.getMeshID("Geometry");
    int forcesID = interface.getDataID("Forces",meshID);

    double data[2] = {5.0,5.0};
    interface.initialize();
    interface.writeVectorData(forcesID,1,data);  //id 1 is hard coded here, maybe, we should get an interface function for this
    interface.advance(0.1);
    interface.finalize();
  }
}


#endif // defined( not PRECICE_NO_MPI )

}} // namespace precice, tests


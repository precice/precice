// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverInterfaceTest.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/config/Configuration.hpp"
#include "geometry/Geometry.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/WrappedVector.h"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerIntegrationTest(precice::tests::SolverInterfaceTest)

namespace precice {
namespace tests {

tarch::logging::Log SolverInterfaceTest:: _log("precice::tests::SolverInterfaceTest");

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
  preciceTrace("run()");
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testConfiguration();
  }
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 1){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
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
      //TODO not working no riemann (benjamin's laptop)
      //testMethod(testBug);
      // testMethod(testNASTINMeshRestart);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  if (Par::getCommunicatorSize() > 2){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 2){
      Par::setGlobalCommunicator(comm);
      testMethod(testThreeSolvers);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  Par::synchronizeProcesses();
  if (Par::getCommunicatorSize() > 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2 , 3;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
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
  preciceTrace1("configureSolverInterface()", configFilename);
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  precice::utils::Events_Clear();
  utils::MasterSlave::reset();
  config::Configuration config;
  utils::configure(config.getXMLTag(), configFilename);
  //validate ( config.isValid() );
  interface._impl->configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceTest:: testConfiguration()
{
  preciceTrace("testConfiguration");
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
  validateEquals(meshContexts[0], static_cast<void*>(NULL));
  validateEquals(meshContexts[1]->mesh->getName(), std::string("ComsolNodes"));
    // Can be replaced by nullptr as soon as we have C++11 available.
  validateEquals(meshContexts[2], static_cast<void*>(NULL));
  validateEquals(comsol->_usedMeshContexts.size(), 1);
}

void SolverInterfaceTest:: testExplicit()
{
  preciceTrace ( "testExplicit" );
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
  preciceTrace("testExplicitWithSubcycling()");
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
  preciceTrace("testExplicitWithDataExchange()");
  using namespace tarch::la;
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  double counter = 0.0;
  using utils::Vector3D;

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
    Vector3D vertex(0.0);
    cplInterface.setMeshVertex(meshOneID, raw(vertex));
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
        Vector3D force(Vector3D(counter) + wrap<3,double>(it.vertexCoords()));
        cplInterface.writeVectorData(forcesID, it.vertexID(), raw(force));
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        i=0;
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          Vector3D vel ( 0.0 );
          int index = indices[i];
          i++;
          cplInterface.readVectorData(velocitiesID, index, raw(vel));
          validate(equals(vel, Vector3D(counter) + wrap<3,double>(it.vertexCoords())));
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
      Vector3D force ( 0.0 );
      cplInterface.readVectorData(forcesID, it.vertexID(), raw(force));
      validate(equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())));
    }
    counter += 1.0;

    while (cplInterface.isCouplingOngoing()){
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3D vel(Vector3D(counter - 1.0) + wrap<3,double>(it.vertexCoords()));
        cplInterface.writeVectorData(velocitiesID, it.vertexID(), raw(vel));
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          Vector3D force(0.0);
          cplInterface.readVectorData(forcesID, it.vertexID(), raw(force));
          validate(equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithDataInitialization()
{
  preciceTrace("testExplicitWithDataInitialization()");
  using namespace tarch::la;
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  using utils::Vector3D;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface cplInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-data-init.xml", cplInterface);
    int meshOneID = cplInterface.getMeshID("MeshOne");
    cplInterface.setMeshVertex(meshOneID, raw(Vector3D(1.0,2.0,3.0)));
    double maxDt = cplInterface.initialize();
    int dataAID = cplInterface.getDataID("DataOne",meshOneID);
    int dataBID = cplInterface.getDataID("DataTwo",meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    validateNumericalEquals(2.0, valueDataB);
    while (cplInterface.isCouplingOngoing()){
      Vector3D valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, raw(valueDataA));
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
    Vector3D pos(0.0);
    cplInterface.setMeshVertex(meshTwoID, raw(pos));
    double maxDt = cplInterface.initialize();
    int dataAID = cplInterface.getDataID("DataOne",meshTwoID);
    int dataBID = cplInterface.getDataID("DataTwo",meshTwoID);
    cplInterface.writeScalarData(dataBID, 0, 2.0);
    //sagen dass daten jetzt geschrieben
    cplInterface.fulfilledAction(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    Vector3D valueDataA;
    cplInterface.readVectorData(dataAID, 0, raw(valueDataA));
    Vector3D expected(1.0, 1.0, 1.0);
    validateWithParams2(equals(valueDataA, expected), valueDataA, expected);
    while (cplInterface.isCouplingOngoing()){
      cplInterface.writeScalarData(dataBID, 0, 2.5);
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, raw(valueDataA));
      validateWithParams2(equals(valueDataA, expected), valueDataA, expected);
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithBlockDataExchange()
{
  preciceTrace("testExplicitWithBlockDataExchange()");
  using namespace tarch::la;
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  double counter = 0.0;
  using utils::Vector3D;

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
    DynamicVector<double> writePositions(size*3);
    DynamicVector<double> getWritePositions(size*3);
    DynamicVector<double> forces(size*3);
    DynamicVector<double> pressures(size);
    DynamicVector<int> writeIDs(size);
    DynamicVector<int> getWriteIDs(size);
    DynamicVector<double> readPositions(size*3);
    DynamicVector<double> getReadPositions(size*3);
    DynamicVector<double> velocities(size*3);
    DynamicVector<double> temperatures(size);
    DynamicVector<double> expectedVelocities(size*3);
    DynamicVector<double> expectedTemperatures(size);
    DynamicVector<int> readIDs(size);
    DynamicVector<int> getReadIDs(size);

    while (cplInterface.isCouplingOngoing()){
      cplInterface.resetMesh(meshOneID);
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        for (int dim=0; dim < 3; dim++){
          writePositions[it.vertexID()*3 + dim] = it.vertexCoords()[dim];
        }
      }
      cplInterface.setMeshVertices(meshOneID, size, raw(writePositions),
                                     raw(writeIDs));
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3D force ( Vector3D(counter) + wrap<3,double>(it.vertexCoords()) );
        for (int dim=0; dim<3; dim++) forces[it.vertexID()*3+dim] = force[dim];
        pressures[it.vertexID()] = counter + it.vertexCoords()[0];
      }
      cplInterface.writeBlockVectorData(forcesID, size, raw(writeIDs), raw(forces));
      cplInterface.writeBlockScalarData(pressuresID, size, raw(writeIDs), raw(pressures));

      cplInterface.getMeshVertices(meshOneID, size, raw(writeIDs),
                                     raw(getWritePositions));
      validateWithParams2(equals(writePositions, getWritePositions),
                          writePositions, getWritePositions);

      cplInterface.getMeshVertexIDsFromPositions(meshOneID, size, raw(writePositions),
                                           raw(getWriteIDs));
      validateWithParams2(equals(writeIDs, getWriteIDs), writeIDs, getWriteIDs);
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
        cplInterface.setMeshVertices(meshOneID, size, raw(readPositions), raw(readIDs));
        cplInterface.mapReadDataTo(meshOneID);
        cplInterface.readBlockVectorData(velocitiesID, size, raw(readIDs),
                                         raw(velocities));
        cplInterface.readBlockScalarData(temperaturesID, size, raw(readIDs),
                                         raw(temperatures));
        validateWithParams2(equals(velocities, expectedVelocities),
                            velocities, expectedVelocities);
        validateWithParams2(equals(temperatures, expectedTemperatures),
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
      Vector3D force ( 0.0 );
      double pressure = 0.0;
      cplInterface.readVectorData(forcesID, it.vertexID(), raw(force));
      cplInterface.readScalarData(pressuresID, it.vertexID(), pressure);
      validate(equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())) );
      validate(equals(pressure, counter + it.vertexCoords()[0]));
    }
    counter += 1.0;

    while ( cplInterface.isCouplingOngoing() ) {
      for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
        Vector3D vel ( Vector3D(counter - 1.0) + wrap<3,double>(it.vertexCoords()) );
        cplInterface.writeVectorData ( velocitiesID, it.vertexID(), raw(vel) );
        double temperature = counter - 1.0 + it.vertexCoords()[0];
        cplInterface.writeScalarData(temperaturesID, it.vertexID(), temperature);
      }
      maxDt = cplInterface.advance ( maxDt );
      if ( cplInterface.isCouplingOngoing() ) {
        for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
          Vector3D force ( 0.0 );
          double pressure = 0.0;
          cplInterface.readVectorData(forcesID, it.vertexID(), raw(force));
          cplInterface.readScalarData(pressuresID, it.vertexID(), pressure);
          validate ( equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())) );
          validate(equals(pressure, counter + it.vertexCoords()[0]));
        }
        counter += 1.0;
      }
    }
    cplInterface.finalize();
  }
}

void SolverInterfaceTest:: testExplicitWithSolverGeometry ()
{
  preciceTrace ( "testExplicitWithSolverGeometry()" );
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
    int i0 = cplInterface.setMeshVertex(meshID, raw(utils::Vector3D(0.0,0.0,0.0)));
    int i1 = cplInterface.setMeshVertex(meshID, raw(utils::Vector3D(1.0,0.0,0.0)));
    int i2 = cplInterface.setMeshVertex(meshID, raw(utils::Vector3D(0.0,1.0,0.0)));
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
   preciceTrace ( "testExplicitWithDisplacingGeometry()" );
   assertion ( utils::Parallel::getCommunicatorSize() > 1 );

   using namespace tarch::la;
   using utils::Vector3D;
   int timesteps = 0;
   double time = 0;

   Vector3D zero (0.0);
   Vector3D one ( zero );
   one(0) += 1.0;
   Vector3D two ( zero );
   two(1) += 1.0;
   Vector3D three ( zero + Vector3D(1.0) );

   if ( utils::Parallel::getProcessRank() == 0 ) { // SolverOne part
      // SolverOne has a local offset to the geometry
      Vector3D localOffset ( 10.0 );
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
         validate(equals(wrap<3,double>(iter.vertexCoords()), zero));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), one));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), two));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), three));
         iter++;
         validate(not (iter != handle.vertices().end()) );

         time += dt;
         timesteps++;
         dt = cplInterface.advance(dt);

         // Add displacements (known in this test case) for validation
         zero += Vector3D(1.0);
         one += Vector3D(2.0);
         two += Vector3D(3.0);
         three += Vector3D(4.0);
      }
      cplInterface.finalize();
   }
   else if ( utils::Parallel::getProcessRank() == 1 ){
      SolverInterface cplInterface ( "SolverTwo", 0, 1 );
      configureSolverInterface (
          _pathToTests + "/explicit-solvergeometry.xml",
          cplInterface );
      int meshID = cplInterface.getMeshID("SolverGeometry");
      int i0 = cplInterface.setMeshVertex(meshID, raw(zero));
      int i1 = cplInterface.setMeshVertex(meshID, raw(one));
      int i2 = cplInterface.setMeshVertex(meshID, raw(two));
      int i3 = cplInterface.setMeshVertex(meshID, raw(three));
//      cplInterface.setMeshEdge ( meshID, i0, i1 );
//      cplInterface.setMeshEdge ( meshID, i1, i3 );
//      cplInterface.setMeshEdge ( meshID, i3, i2 );
//      cplInterface.setMeshEdge ( meshID, i2, i0 );

      double dt = cplInterface.initialize();
      int displacementsID = cplInterface.getDataID("Displacements", meshID);

      while (cplInterface.isCouplingOngoing()){
         MeshHandle handle = cplInterface.getMeshHandle("SolverGeometry");
         VertexIterator iter = handle.vertices().begin();
         validate(equals(wrap<3,double>(iter.vertexCoords()), zero));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), one));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), two));
         iter++;
         validate(equals(wrap<3,double>(iter.vertexCoords()), three));
         iter++;
         validate(not (iter != handle.vertices().end()));

         // Add displacements
         cplInterface.writeVectorData(displacementsID, i0, raw(Vector3D(1.0)));
         cplInterface.writeVectorData(displacementsID, i1, raw(Vector3D(2.0)));
         cplInterface.writeVectorData(displacementsID, i2, raw(Vector3D(3.0)));
         cplInterface.writeVectorData(displacementsID, i3, raw(Vector3D(4.0)));

         // modify coordinates by displacements
         zero += Vector3D(1.0);
         one += Vector3D(2.0);
         two += Vector3D(3.0);
         three += Vector3D(4.0);

         time += dt;
         dt = cplInterface.advance(dt);
         timesteps++;
      }
      cplInterface.finalize();
   }
}

void SolverInterfaceTest:: testExplicitWithDataScaling()
{
  preciceTrace ( "testExplicitWithDataScaling" );
  assertion ( utils::Parallel::getCommunicatorSize() == 2 );
  using namespace tarch::la;
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
        utils::Vector2D data((double)i);
        cplInterface.writeVectorData ( velocitiesID, i, raw(data) );
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
        utils::Vector2D readData;
        cplInterface.readVectorData ( velocitiesID, i, raw(readData) );
        utils::Vector2D expectedData ( (double)i * 10.0 );
        validate ( equals(readData, expectedData, 5e-13) );
      }
      dt = cplInterface.advance(dt);
    }
    cplInterface.finalize();
  }
}


void SolverInterfaceTest:: testExplicitWithCheckpointingStatMapping()
{
  preciceTrace("testExplicitWithCheckpointingStatMapping()");
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  using namespace tarch::la;
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
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->fromData;
      validateEquals(forcesID, dataContext->fromData->getID());
      utils::DynVector& values = data->values();
      assign(values) = 1.0;
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
      validate(dataContext != NULL);
      mesh::PtrData fromData = dataContext->fromData;
      mesh::PtrData toData = dataContext->toData;
      // no mapping, so fromData and toData should be the same data
      validateEquals(fromData->getID(), toData->getID());
      validateEquals(dataID, fromData->getID());
      utils::DynVector& values = fromData->values();
      assign(values) = 2.0;
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
    utils::Vector2D integral(0.0);
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
    utils::Vector2D integral(0.0);
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
  preciceTrace("testImplicit");
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
  preciceTrace("testImplicitWithCheckpointingMappingStat()");
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
      validate(dataContext != NULL);
      mesh::PtrData localForces = dataContext->fromData;
      validateEquals(forcesID, localForces->getID());
      utils::DynVector& forceValues = localForces->values();
      assign(forceValues) = 1.0;
      iterationCount++;
      maxDt = couplingInterface.advance(maxDt);

      dataContext = solverOne->_dataContexts[velocitiesID];
      validate(dataContext != NULL);
      mesh::PtrData localVelocities = dataContext->toData;
      validateEquals(velocitiesID, localVelocities->getID());
      utils::DynVector& velocityValues = localVelocities->values();
      utils::Vector2D integral;
      sumSubvectors(velocityValues, integral);
      validate(equals(integral, utils::Vector2D(8.0, 8.0)));

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
      validate(dataContext != NULL);
      mesh::PtrData velocities = dataContext->fromData;
      validateEquals(velocitiesID, velocities->getID());
      utils::DynVector& values = velocities->values();
      assign(values) = 2.0;

      maxDt = couplingInterface.advance(maxDt);

      dataContext = solverTwo->_dataContexts[forcesID];
      validate(dataContext != NULL);
      mesh::PtrData forces = dataContext->toData;
      //validateEquals(velocitiesID, localVelocities->getID());
      utils::DynVector& forceValues = forces->values();
      utils::Vector2D integral;
      sumSubvectors(forceValues, integral);
      validate(equals(integral, utils::Vector2D(4.0, 4.0)));

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
    validate(dataContext != NULL);
    mesh::PtrData localVelocities = dataContext->toData;
    validateEquals(velocitiesID, localVelocities->getID());
    utils::DynVector& velocityValues = localVelocities->values();
    utils::Vector2D integral;
    sumSubvectors(velocityValues, integral);
    validate(equals(integral, utils::Vector2D(8.0, 8.0)));

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
    validate(dataContext != NULL);
    mesh::PtrData forces = dataContext->toData;
    validateEquals(forcesID, forces->getID());
    utils::DynVector& forceValues = forces->values();
    utils::Vector2D integral;
    sumSubvectors(forceValues, integral);
    validate(equals(integral, utils::Vector2D(4.0, 4.0)));

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
  preciceTrace2("runSolver()", solverName, configurationFileName);
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
  preciceTrace("testStationaryMappingWithSolverMesh()");
  std::string config2D = _pathToTests + "mapping-without-geo-2D.xml";
  std::string config3D = _pathToTests + "mapping-without-geo-3D.xml";
  int rank = utils::Parallel::getProcessRank();
  assertion1((rank == 0) || (rank == 1), rank);
  std::string solverName = rank == 0 ? "SolverA" : "SolverB";
  std::string meshForcesA = "MeshForcesA";
  std::string meshDisplA = "MeshDisplacementsA";
  std::string meshForcesB = "MeshForcesB";
  std::string meshDisplB = "MeshDisplacementsB";
  std::string dataForces = constants::dataForces();
  std::string dataDispl = constants::dataDisplacements();
  using tarch::la::raw;
  using tarch::la::equals;

  for (int dim=2; dim < 3; dim++){
    preciceDebug("Running " << dim << "D test");
    SolverInterface interface(solverName, 0, 1);
    if (dim == 2){
      configureSolverInterface(config2D, interface);
    }
    else {
      configureSolverInterface(config3D, interface);
    }
    validateEquals(interface.getDimensions(), dim);

    std::vector<utils::DynVector> positions;
    utils::DynVector position(dim);
    if (dim == 2){
      assignList(position) = 0.0, 0.0;
      positions.push_back(position);
      assignList(position) = 1.0, 0.0;
      positions.push_back(position);
      assignList(position) = 1.0, 1.0;
      positions.push_back(position);
      assignList(position) = 0.0, 1.0;
      positions.push_back(position);
    }
    else {
      assignList(position) = 0.0, 0.0, 0.0;
      positions.push_back(position);
      assignList(position) = 1.0, 0.0, 0.0;
      positions.push_back(position);
      assignList(position) = 1.0, 1.0, 0.0;
      positions.push_back(position);
      assignList(position) = 0.0, 1.0, 1.0;
      positions.push_back(position);
      assignList(position) = 0.0, 0.0, 1.0;
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
        interface.setMeshVertex(meshForcesID, raw(positions[i] + 0.1));
        interface.setMeshVertex(meshDisplID, raw(positions[i] + 0.6));
      }
      double maxDt = interface.initialize();

      validate(interface.isWriteDataRequired(maxDt));
      validate(not interface.isReadDataAvailable());
      utils::DynVector force(dim, 1.0);
      utils::DynVector displ(dim, 0.0);
      for (size_t i=0; i < size; i++){
        interface.writeVectorData(dataForcesID, i, raw(force));
      }
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      interface.mapReadDataTo(meshDisplID);
      //precicePrint("1: mapped data: " << interface._impl->_accessor->dataContext(dataDisplID).data->values());
      force += 1.0;
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, raw(displ));
        validateNumericalEquals(displ[0], positions[i][0] + 0.1);
        interface.writeVectorData(dataForcesID, i, raw(force));
      }
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      interface.mapReadDataTo(meshDisplID);
      //precicePrint("2: mapped data: " << interface._impl->_accessor->dataContext(dataDisplID).data->values());
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, raw(displ));
        validateNumericalEquals(displ[0], 2.0*(positions[i][0] + 0.1));
      }
      interface.finalize();
    }
    else {
      assertion1(rank == 1, rank);
      int meshForcesID = interface.getMeshID(meshForcesB);
      int meshDisplID = interface.getMeshID(meshDisplB);
      int dataForcesID = interface.getDataID(dataForces, meshForcesID);
      int dataDisplID = interface.getDataID(dataDispl, meshDisplID);

      // Set solver mesh positions provided to SolverA for data mapping
      for (size_t i=0; i < size; i++){
        interface.setMeshVertex(meshForcesID, raw(positions[i]));
        interface.setMeshVertex(meshDisplID, raw(positions[i] + 0.5));
      }
      double maxDt = interface.initialize();

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      utils::DynVector force(dim, 0.0);
      utils::DynVector totalForce(dim, 0.0);
      utils::DynVector displ(dim, 0.0);
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataForcesID, i, raw(force));
        totalForce += force;
        assign(displ) = positions[i][0];
        interface.writeVectorData(dataDisplID, i, raw(displ));
      }
      utils::DynVector expected(dim, (double)size);
      validateWithParams2(equals(totalForce,expected), totalForce, expected);
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      assign(totalForce) = 0.0;
      for (size_t i=0; i < positions.size(); i++){
        interface.readVectorData(dataForcesID, i, raw(force));
        totalForce += force;
        assign(displ) = 2.0 * positions[i][0];
        interface.writeVectorData(dataDisplID, i, raw(displ));
      }
      assign(expected) = 2.0 * (double)size;
      validateWithParams1(equals(totalForce,expected), totalForce);
      maxDt = interface.advance(maxDt);

      validate(interface.isWriteDataRequired(maxDt));
      validate(interface.isReadDataAvailable());
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, raw(force));
      }
      interface.finalize();
    }
  }
}

void SolverInterfaceTest:: testDistributedCommunications()
{
  preciceTrace("testDistributedCommunications()");

  assertion(utils::Parallel::getCommunicatorSize() == 4);

  std::vector<std::string> fileNames({
      "point-to-point-sockets.xml",
      "point-to-point-mpi.xml",
      "gather-scatter-mpi.xml"});

  for (auto fileName : fileNames) {
    mesh::Mesh::resetGeometryIDsGlobally();

    std::string solverName;
    int rank, size;
    std::string meshName;
    int i1,i2; //indices for data and positions

    std::vector<utils::DynVector> positions;
    std::vector<utils::DynVector> data;
    std::vector<utils::DynVector> expectedData;

    utils::DynVector position(3);
    utils::DynVector datum(3);

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
      int vertexID = precice.setMeshVertex(meshID, raw(positions[i]));
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    if (utils::Parallel::getProcessRank() <= 1){ //Fluid
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.writeVectorData(forcesID, vertexIDs[i], raw(data[i+i1]));
      }
    }
    else if (utils::Parallel::getProcessRank() >= 2){ //Structure
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.readVectorData(forcesID, vertexIDs[i], raw(data[i]));
        data[i] = data[i]*2 + 1.0;
        precice.writeVectorData(velocID, vertexIDs[i], raw(data[i]));
      }
    }

    precice.advance(1.0);

    if (utils::Parallel::getProcessRank() <= 1){ //Fluid
      for( size_t i=0; i<vertexIDs.size(); i++){
        precice.readVectorData(velocID, vertexIDs[i], raw(data[i+i1]));
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
  preciceTrace("testBug()");
  typedef utils::Vector3D Vector3D;
  using namespace tarch::la;
  std::string config = _pathToTests + "bug.xml";

  int slices = 5;
  std::vector<Vector3D> coords;
  for (int i=0; i < slices; i++){
    double z = (double)i * 1.0;
    coords += Vector3D( 1.0,  0.0, z),
              Vector3D( 0.0,  1.0, z),
              Vector3D(-1.0,  0.0, z),
              Vector3D( 0.0, -1.0, z);
  }

  int rank = utils::Parallel::getProcessRank();
  assertion1((rank == 0) || (rank == 1), rank);
  std::string solverName = rank == 0 ? "Flite" : "Calculix";
  if (solverName == std::string("Flite")){
    SolverInterface precice("Flite", 0, 1);
    configureSolverInterface(config, precice);
    int meshID = precice.getMeshID("FliteNodes");
    int forcesID = precice.getDataID(precice::constants::dataForces(), meshID);
    int displacementsID = precice.getDataID(precice::constants::dataDisplacements(), meshID);
    int oldDisplacementsID = precice.getDataID("OldDisplacements", meshID);
    validateEquals(precice.getDimensions(), 3);
    foreach(Vector3D& coord, coords){
      precice.setMeshVertex(meshID, raw(coord));
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
    assertion1(solverName == std::string("Calculix"), solverName);
    SolverInterface precice("Calculix", 0, 1);
    configureSolverInterface(config, precice);
    int meshID = precice.getMeshID("CalculixNodes");
    foreach (Vector3D& coord, coords){
      precice.setMeshVertex(meshID, raw(coord));
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
  preciceTrace("testThreeSolvers()");
  std::string configFilename(_pathToTests + "three-solver-explicit-explicit.xml");
  std::vector<int> expectedCallsOfAdvance;
  expectedCallsOfAdvance += 10, 10, 10;
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-implicit-implicit.xml";
  expectedCallsOfAdvance.clear();
  expectedCallsOfAdvance += 30, 30, 20;
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-implicit-explicit.xml";
  expectedCallsOfAdvance.clear();
  expectedCallsOfAdvance += 30, 30, 10;
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-explicit-implicit.xml";
  expectedCallsOfAdvance.clear();
  expectedCallsOfAdvance += 30, 10, 30;
  runThreeSolvers(configFilename, expectedCallsOfAdvance);

  configFilename = _pathToTests + "three-solver-parallel.xml";
  expectedCallsOfAdvance.clear();
  expectedCallsOfAdvance += 30, 30, 10;
  runThreeSolvers(configFilename, expectedCallsOfAdvance);
}

void SolverInterfaceTest:: runThreeSolvers
(
  const std::string&      configFilename,
  const std::vector<int>& expectedCallsOfAdvance )
{
  preciceTrace2("runThreeSolvers", configFilename, expectedCallsOfAdvance);

  typedef utils::Vector2D Vector2D;
  int rank = utils::Parallel::getProcessRank();
  assertion1((rank == 0) || (rank == 1) || (rank == 2), rank);

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
    precice.setMeshVertex(meshID, raw(utils::Vector2D(0.0, 0.0)));
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
    assertion1(solverName == std::string("SolverThree"), solverName);
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
  preciceTrace("testMultiCoupling()");
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  mesh::Mesh::resetGeometryIDsGlobally();

  std::vector<utils::DynVector> positions;
  utils::DynVector position(2);
  assignList(position) = 0.0, 0.0;
  positions.push_back(position);
  assignList(position) = 1.0, 0.0;
  positions.push_back(position);
  assignList(position) = 1.0, 1.0;
  positions.push_back(position);
  assignList(position) = 0.0, 1.0;
  positions.push_back(position);

  std::vector<utils::DynVector> datas;
  utils::DynVector data(2);
  assignList(data) = 1.0, 1.0;
  datas.push_back(data);
  assignList(data) = 2.0, 2.0;
  datas.push_back(position);
  assignList(data) = 3.0, 3.0;
  datas.push_back(data);
  assignList(data) = 4.0, 5.0;
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
      vertexID = precice.setMeshVertex(meshID, raw(positions[i]));
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i=0; i < 4; i++){
      precice.writeVectorData(dataWriteID, vertexIDs[i], raw(datas[i]));
    }

    if (precice.isActionRequired(writeIterCheckpoint)){
      precice.fulfilledAction(writeIterCheckpoint);
    }
    precice.advance(0.0001);
    if (precice.isActionRequired(readIterCheckpoint)){
      precice.fulfilledAction(readIterCheckpoint);
    }

    for (size_t i=0; i < 4; i++){
      precice.readVectorData(dataReadID, vertexIDs[i], raw(datas[i]));
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
      vertexID = precice.setMeshVertex(meshID1, raw(positions[i]));
      vertexIDs1.push_back(vertexID);
    }
    std::vector<int> vertexIDs2;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID2, raw(positions[i]));
      vertexIDs2.push_back(vertexID);
    }
    std::vector<int> vertexIDs3;
    for (size_t i=0; i < 4; i++){
      vertexID = precice.setMeshVertex(meshID3, raw(positions[i]));
      vertexIDs3.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i=0; i < 4; i++){
      precice.writeVectorData(dataWriteID1, vertexIDs1[i], raw(datas[i]));
      precice.writeVectorData(dataWriteID2, vertexIDs2[i], raw(datas[i]));
      precice.writeVectorData(dataWriteID3, vertexIDs3[i], raw(datas[i]));
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
  preciceTrace("testNASTINMeshRestart()");
  assertion(utils::Parallel::getCommunicatorSize() == 2);

  std::vector<std::string> restartFiles;
  restartFiles.push_back("precice_checkpoint_NASTIN_NASTIN_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_NASTIN_SOLIDZ_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_NASTIN_simstate.txt");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_cplscheme.txt");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_SOLIDZ_Mesh.wrl");
  restartFiles.push_back("precice_checkpoint_SOLIDZ_simstate.txt");

  foreach(std::string& restartFile, restartFiles){
    std::ifstream  src((_pathToTests + restartFile).c_str(), std::ifstream::in);
    std::ofstream  dst(restartFile.c_str(), std::ifstream::out);
    dst << src.rdbuf();
  }

  std::string readSimCheckpoint(constants::actionReadSimulationCheckpoint());
  std::string writeItCheckpoint(constants::actionWriteIterationCheckpoint());

  mesh::Mesh::resetGeometryIDsGlobally();
  int meshSize = 27;

  double positions[meshSize*2];

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
  precice.setMeshVertices(meshID, meshSize, positions, vertexIDs);
  precice.initialize();

  if (precice.isActionRequired(writeItCheckpoint)){
    precice.fulfilledAction(writeItCheckpoint);
  }
  precice.finalize();
}

#endif // defined( not PRECICE_NO_MPI )

}} // namespace precice, tests

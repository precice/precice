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
      testMethod(testExplicitWithBlockDataExchange);
      testMethod(testExplicitWithSolverGeometry);
      testMethod(testExplicitWithDisplacingGeometry);
      testMethod(testExplicitWithDataScaling);
#     ifndef PRECICE_NO_SPIRIT2
      testMethod(testExplicitWithCheckpointingIncMapping);
      testMethod(testExplicitWithCheckpointingStatMapping);
#     endif // not PRECICE_NO_SPIRIT2
      testMethod(testImplicit);
#     ifndef PRECICE_NO_SPIRIT2
      testMethod(testImplicitWithCheckpointingMappingInc);
      testMethod(testImplicitWithCheckpointingMappingStat);
#     endif // not PRECICE_NO_SPIRIT2
      testMethod(testStationaryMappingWithSolverMesh);
      testMethod(testBug);
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
  validateEquals ( meshContexts.size(), 2 );
  validateEquals ( peano->_usedMeshContexts.size(), 2 );

  validateEquals ( meshContexts[0]->mesh->getName(), std::string("Sphere") );
  validateEquals ( meshContexts[1]->mesh->getName(), std::string("Boundingbox") );

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
  validateEquals(meshContexts.size(), 2);
  validateEquals(meshContexts[0]->mesh->getName(), std::string("Sphere"));
  validateEquals(meshContexts[1], NULL);
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
    double maxDt = cplInterface.initialize();
    int forcesID = cplInterface.getDataID("Forces");
    int velocitiesID = cplInterface.getDataID("Velocities");
    int meshID = cplInterface.getMeshID("Test-Square");
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    while (cplInterface.isCouplingOngoing()){
      cplInterface.resetWritePositions(meshID);
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        int index = cplInterface.setWritePosition(meshID, it.vertexCoords());
        validateEquals(index, it.vertexID());
      }
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3D force(Vector3D(counter) + wrap<3,double>(it.vertexCoords()));
        cplInterface.writeVectorData(forcesID, it.vertexID(), raw(force));
      }
      maxDt = cplInterface.advance(maxDt);
      if (cplInterface.isCouplingOngoing()){
        for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
          Vector3D vel ( 0.0 );
          int index = cplInterface.setReadPosition(meshID, it.vertexCoords());
          validateEquals(index, 0);
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
    int forcesID = cplInterface.getDataID("Forces");
    int velocitiesID = cplInterface.getDataID("Velocities");
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
    int forcesID = cplInterface.getDataID("Forces");
    int velocitiesID = cplInterface.getDataID("Velocities");
    int meshID = cplInterface.getMeshID("Test-Square");
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    int size = vertices.size();
    DynamicVector<double> writePositions(size*3);
    DynamicVector<double> getWritePositions(size*3);
    DynamicVector<double> forces(size*3);
    DynamicVector<int> writeIDs(size);
    DynamicVector<int> getWriteIDs(size);
    DynamicVector<double> readPositions(size*3);
    DynamicVector<double> getReadPositions(size*3);
    DynamicVector<double> velocities(size*3);
    DynamicVector<double> expectedVelocities(size*3);
    DynamicVector<int> readIDs(size);
    DynamicVector<int> getReadIDs(size);

    while (cplInterface.isCouplingOngoing()){
      cplInterface.resetWritePositions(meshID);
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        for (int dim=0; dim < 3; dim++){
          writePositions[it.vertexID()*3 + dim] = it.vertexCoords()[dim];
        }
      }
      cplInterface.setWritePositions(meshID, size, raw(writePositions),
                                     raw(writeIDs));
      for (VertexIterator it = vertices.begin(); it != vertices.end(); it++){
        Vector3D force ( Vector3D(counter) + wrap<3,double>(it.vertexCoords()) );
        for (int dim=0; dim<3; dim++) forces[it.vertexID()*3+dim] = force[dim];
      }
      cplInterface.writeBlockVectorData(forcesID, size, raw(writeIDs), raw(forces));

      cplInterface.getWritePositions(meshID, size, raw(writeIDs),
                                     raw(getWritePositions));
      validateWithParams2(equals(writePositions, getWritePositions),
                          writePositions, getWritePositions);

      cplInterface.getWriteIDsFromPositions(meshID, size, raw(writePositions),
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

        }
        cplInterface.resetReadPositions(meshID);
        cplInterface.setReadPositions(meshID, size, raw(readPositions), raw(readIDs));
        cplInterface.mapReadData(meshID);
        cplInterface.readBlockVectorData(velocitiesID, size, raw(readIDs),
                                         raw(velocities));
        validateWithParams2(equals(velocities, expectedVelocities),
                            velocities, expectedVelocities);

        cplInterface.getReadPositions(meshID, size, raw(readIDs),
                                      raw(getReadPositions));
        validateWithParams2(equals(readPositions, getReadPositions),
                            readPositions, getReadPositions);

        cplInterface.getReadIDsFromPositions(meshID, size, raw(readPositions),
                                             raw(getReadIDs));
        validateWithParams2(equals(readIDs, getReadIDs), readIDs, getReadIDs);
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
    int forcesID = cplInterface.getDataID ( "Forces" );
    int velocitiesID = cplInterface.getDataID ( "Velocities" );
    VertexHandle vertices = cplInterface.getMeshHandle("Test-Square").vertices();
    // SolverTwo does not start the coupled simulation and has, hence,
    // already received the first data to be validated.
    for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
      Vector3D force ( 0.0 );
      cplInterface.readVectorData ( forcesID, it.vertexID(), raw(force) );
      validate ( equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())) );
    }
    counter += 1.0;

    while ( cplInterface.isCouplingOngoing() ) {
      for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
        Vector3D vel ( Vector3D(counter - 1.0) + wrap<3,double>(it.vertexCoords()) );
        cplInterface.writeVectorData ( velocitiesID, it.vertexID(), raw(vel) );
      }
      maxDt = cplInterface.advance ( maxDt );
      if ( cplInterface.isCouplingOngoing() ) {
        for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
          Vector3D force ( 0.0 );
          cplInterface.readVectorData ( forcesID, it.vertexID(), raw(force) );
          validate ( equals(force, Vector3D(counter) + wrap<3,double>(it.vertexCoords())) );
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

  int timesteps;
  double time;

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
    int i0 = cplInterface.setMeshVertex ( meshID, raw(utils::Vector3D(1.0)) );
    int i1 = cplInterface.setMeshVertex ( meshID, raw(utils::Vector3D(2.0)) );
    cplInterface.setMeshEdge ( meshID, i0, i1 );
    double dt = cplInterface.initialize();
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
   int timesteps;
   double time;

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
      while ( cplInterface.isCouplingOngoing() ){
         MeshHandle handle = cplInterface.getMeshHandle ( "SolverGeometry" );
         VertexIterator iter = handle.vertices().begin();
         validate ( equals(wrap<3,double>(iter.vertexCoords()), zero) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), one) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), two) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), three) );
         iter++;
         validate ( not (iter != handle.vertices().end()) );

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
      int meshID = cplInterface.getMeshID ( "SolverGeometry" );
      int i0 = cplInterface.setMeshVertex ( meshID, raw(zero) );
      int i1 = cplInterface.setMeshVertex ( meshID, raw(one) );
      int i2 = cplInterface.setMeshVertex ( meshID, raw(two) );
      int i3 = cplInterface.setMeshVertex ( meshID, raw(three) );
      cplInterface.setMeshEdge ( meshID, i0, i1 );
      cplInterface.setMeshEdge ( meshID, i1, i3 );
      cplInterface.setMeshEdge ( meshID, i3, i2 );
      cplInterface.setMeshEdge ( meshID, i2, i0 );

      double dt = cplInterface.initialize();
      int displacementsID = cplInterface.getDataID ( "Displacements" );

      while ( cplInterface.isCouplingOngoing() ){
         MeshHandle handle = cplInterface.getMeshHandle ( "SolverGeometry" );
         VertexIterator iter = handle.vertices().begin();
         validate ( equals(wrap<3,double>(iter.vertexCoords()), zero) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), one) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), two) );
         iter++;
         validate ( equals(wrap<3,double>(iter.vertexCoords()), three) );
         iter++;
         validate ( not (iter != handle.vertices().end()) );

         // Add displacements
         cplInterface.writeVectorData ( displacementsID, i0, raw(Vector3D(1.0)) );
         cplInterface.writeVectorData ( displacementsID, i1, raw(Vector3D(2.0)) );
         cplInterface.writeVectorData ( displacementsID, i2, raw(Vector3D(3.0)) );
         cplInterface.writeVectorData ( displacementsID, i3, raw(Vector3D(4.0)) );

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

//void SolverInterfaceTest:: testCoupledSimulationWithPreCoupledTimesteps ()
//{
//   preciceDebug ( "testCoupledSimulationWithPreCoupledTimesteps()",
//                   "Entering" );
//
//   int timesteps = 0;
//   double time = 0.0;
//
//   if ( utils::Parallel::getProcessRank() == 0 ) {
//      runSolver ( "SolverOne",   preciceDebug ( "testTrivialCoupledSimulationWithDataExchange", "Entering" );
//                  "/SolverInterfaceTest-precomputations-config.xml",
//                  timesteps, time );
//      validateEquals ( timesteps, 20,
//                       "testCoupledSimulationWithPreCoupledTimesteps" );
//      validateNumericalEquals ( time, 0.2,
//                                "testCoupledSimulationWithPreCoupledTimesteps" );
//   }
//   else if ( utils::Parallel::getProcessRank() == 1 ) {
//      runSolver ( "SolverTwo",
//                  "/SolverInterfaceTest-precomputations-config.xml",
//                  timesteps, time );
//      validateEquals ( timesteps, 10,
//                       "testCoupledSimulationWithPreCoupledTimesteps" );
//      validateNumericalEquals ( time, 0.1,
//                                "testCoupledSimulationWithPreCoupledTimesteps" );
//   }
//
//   preciceDebug ( "testCoupledSimulationWithPreCoupledTimesteps()",
//                   "Leaving" );
//}

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
    int velocitiesID = cplInterface.getDataID ( "Velocities" );
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
    int velocitiesID = cplInterface.getDataID ( "Velocities" );
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

void SolverInterfaceTest:: testExplicitWithCheckpointingIncMapping()
{
  preciceTrace("testExplicitWithCheckpointingIncMapping()");
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::Mesh::resetGeometryIDsGlobally();
  using namespace tarch::la;
  using namespace precice::constants;
  int timesteps = 0;
  double time = 0.0;

  if (utils::Parallel::getProcessRank() == 0){
    SolverInterface couplingInterface("SolverOne", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-writecheckpoint-inc.xml",
                             couplingInterface);
    validateEquals ( couplingInterface.getDimensions(), 2 );
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals ( solverOne->getName(), "SolverOne" );
    double dt = couplingInterface.initialize();
    int dataID = couplingInterface.getDataID("Forces");
    validateEquals(solverOne->_meshContexts.size(), 1);
    while (couplingInterface.isCouplingOngoing()){
      validate(solverOne->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverOne->_dataContexts[dataID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->data;
      validateEquals(dataID, data->getID());
      utils::DynVector& values = data->values();
      assign(values) = 1.0;
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
  else if (utils::Parallel::getProcessRank() == 1){
    SolverInterface couplingInterface("SolverTwo", 0, 1);
    configureSolverInterface(_pathToTests + "/explicit-writecheckpoint-inc.xml",
                             couplingInterface );
    validateEquals(couplingInterface.getDimensions(), 2);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    double dt = couplingInterface.initialize();
    int dataID = couplingInterface.getDataID("Velocities");
    validateEquals(solverTwo->_meshContexts.size(), 1);
    while (couplingInterface.isCouplingOngoing()){
      validate(solverTwo->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverTwo->_dataContexts[dataID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->data;
      validateEquals(dataID, data->getID());
      utils::DynVector& values = data->values();
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
    configureSolverInterface(_pathToTests + "/explicit-readcheckpoint-inc.xml",
                             couplingInterface);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    double dt = couplingInterface.initialize();
    int dataID = couplingInterface.getDataID("Forces");
    validateEquals(solverOne->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverOne->_meshContexts[0]->mesh;
    utils::Vector2D integral(0.0);
    couplingInterface.integrateVectorData(dataID, raw(integral));
    validate(tarch::la::equals(integral, utils::Vector2D(40.0)));
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
    configureSolverInterface(_pathToTests + "/explicit-readcheckpoint-inc.xml",
                             couplingInterface);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    double dt = couplingInterface.initialize();
    int dataID = couplingInterface.getDataID("Velocities");
    validateEquals(solverTwo->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverTwo->_meshContexts[0]->mesh;
    utils::Vector2D integral(0.0);
    couplingInterface.integrateVectorData(dataID, raw(integral));
    validate(tarch::la::equals(integral, utils::Vector2D(80.0)));
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

    int meshID = couplingInterface.getMeshID("Test-Square");
    double pos[2];
    // Set write positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setWritePosition(meshID, pos);
    // Set read positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setReadPosition(meshID, pos);

    double dt = couplingInterface.initialize();
    int forcesID = couplingInterface.getDataID("Forces");
    int velocitiesID = couplingInterface.getDataID("Velocities");
    validateEquals(solverOne->_meshContexts.size(), 1);

    while (couplingInterface.isCouplingOngoing()){
      validate(solverOne->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverOne->_dataContexts[forcesID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->localData;
      validateEquals(forcesID, dataContext->data->getID());
      utils::DynVector& values = data->values();
      assign(values) = 1.0;
      time += dt;
      dt = couplingInterface.advance(dt);
      couplingInterface.mapReadData(meshID);
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
    int dataID = couplingInterface.getDataID("Velocities");
    validateEquals(solverTwo->_meshContexts.size(), 1);
    while (couplingInterface.isCouplingOngoing()){
      validate(solverTwo->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverTwo->_dataContexts[dataID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->data;
      validateEquals(dataID, data->getID());
      utils::DynVector& values = data->values();
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
    double dt = couplingInterface.initialize();
    int forcesID = couplingInterface.getDataID("Forces");
    validateEquals(solverOne->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverOne->_meshContexts[0]->mesh;
    utils::Vector2D integral(0.0);
    couplingInterface.integrateVectorData(forcesID, raw(integral));
    validateWithMessage(tarch::la::equals(integral, utils::Vector2D(4.0)),
                        integral);
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
    int dataID = couplingInterface.getDataID("Velocities");
    validateEquals(solverTwo->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverTwo->_meshContexts[0]->mesh;
    utils::Vector2D integral(0.0);
    couplingInterface.integrateVectorData(dataID, raw(integral));
    validate(tarch::la::equals(integral, utils::Vector2D(8.0)));
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

void SolverInterfaceTest:: testImplicitWithSubcycling()
{
//   preciceTrace ( "testImplicitWithSubcycling()" );
//   assertion ( utils::Parallel::getCommunicatorSize() > 1  );
//
//   int timesteps = 0;
//   double time = 0.0;
//
//   if ( utils::Parallel::getProcessRank() == 0 ) {
//      SolverInterface couplingInterface ( "SolverOne" );
//      couplingInterface.configure (
//         _pathToTests +
//         "/SolverInterfaceTest-implicit-subcycling-config.xml" );
//      double dt = couplingInterface.initialize ();
//      int dataID = couplingInterface.getDataID ( "Forces" );
//      validateEquals ( couplingInterface._meshContexts.size(), 1, "testCheckpoint" );
//      while ( couplingInterface.isCouplingOngoing() ) {
//         mesh::PtrMesh mesh = couplingInterface._meshContexts[0].mesh;
//         foreach ( mesh::Vertex & vertex, mesh->vertices() ) {
//            vertex.setProperty ( dataID, Vector(1.0) );
//         }
//         time += dt;
//         dt = couplingInterface.advance ( dt );
//         timesteps++;
//         if ( couplingInterface.isActionRequired(
//              Constants::NAME_WRITE_SIMULATION_CHECKPOINT) )
//         {
//            validateEquals ( timesteps, 10, "testExplicitWithCheckpointingIncMapping" );
//            couplingInterface.fulfilledAction (
//               Constants::NAME_WRITE_SIMULATION_CHECKPOINT );
//         }
//      }
//      couplingInterface.finalize ();
//
//      validateEquals ( timesteps, 10, "testExplicitWithCheckpointingIncMapping" );
//      validateEquals ( time, 10.0, "testExplicitWithCheckpointingIncMapping" );
//   }
//   else if ( utils::Parallel::getProcessRank() == 1 ) {
//      SolverInterface couplingInterface ( "SolverTwo" );
//      couplingInterface.configure (
//         _pathToTests +
//         "/SolverInterfaceTest-writecheckpoint-config.xml" );
//      double dt = couplingInterface.initialize ();
//      int dataID = couplingInterface.getDataID ( "Velocities" );
//      validateEquals ( couplingInterface._meshContexts.size(), 1, "testCheckpoint" );
//      while ( couplingInterface.isCouplingOngoing() ) {
//         mesh::PtrMesh mesh = couplingInterface._meshContexts[0].mesh;
//         foreach ( mesh::Vertex & vertex, mesh->vertices() ) {
//            vertex.setProperty ( dataID, Vector(2.0) );
//         }
//         time += dt;
//         dt = couplingInterface.advance ( dt );
//         timesteps++;
//         if ( couplingInterface.isActionRequired(
//              Constants::NAME_WRITE_SIMULATION_CHECKPOINT) )
//         {
//            validateEquals ( timesteps, 10, "testExplicitWithCheckpointingIncMapping" );
//            couplingInterface.fulfilledAction (
//               Constants::NAME_WRITE_SIMULATION_CHECKPOINT );
//         }
//      }
//      couplingInterface.finalize ();
//
//      validateEquals ( timesteps, 10, "testExplicitWithCheckpointingIncMapping" );
//      validateEquals ( time, 10.0, "testExplicitWithCheckpointingIncMapping" );
//   }
}

void SolverInterfaceTest:: testImplicitWithCheckpointingMappingInc()
{
  preciceTrace("testImplicitWithCheckpointingMappingInc()");
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
    configureSolverInterface(_pathToTests + "/implicit-writecheckpoint-inc.xml", couplingInterface);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    int dataID = couplingInterface.getDataID("Forces");
    double maxDt = couplingInterface.initialize();

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      validate(solverOne->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverOne->_dataContexts[dataID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->data;
      validateEquals(dataID, data->getID());
      utils::DynVector& values = data->values();
      assign(values) = 1.0;
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
    configureSolverInterface(_pathToTests + "/implicit-writecheckpoint-inc.xml", couplingInterface);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    int dataID = couplingInterface.getDataID("Velocities");
    double maxDt = couplingInterface.initialize();

    while (couplingInterface.isCouplingOngoing()){
      if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())){
        couplingInterface.fulfilledAction(actionWriteIterationCheckpoint());
        checkpoint = state;
        iterationCount = 1;
      }
      validate(solverTwo->_dataContexts.size() > 0);
      impl::DataContext* dataContext = solverTwo->_dataContexts[dataID];
      validate(dataContext != NULL);
      mesh::PtrData data = dataContext->data;
      validateEquals(dataID, data->getID());
      utils::DynVector& values = data->values();
      assign(values) = 2.0;
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
    configureSolverInterface(_pathToTests + "/implicit-readcheckpoint-inc.xml", couplingInterface);
    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());
    int dataID = couplingInterface.getDataID("Forces");
    validateEquals(solverOne->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverOne->_meshContexts[0]->mesh;
    utils::Vector3D integral(0.0);
    couplingInterface.integrateVectorData(dataID, raw(integral));
    validateWithMessage(tarch::la::equals(integral, utils::Vector3D(8.0)), integral);

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
    configureSolverInterface(_pathToTests + "/implicit-readcheckpoint-inc.xml", couplingInterface);
    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());
    int dataID = couplingInterface.getDataID("Velocities");
    validateEquals(solverTwo->_meshContexts.size(), 1);
    mesh::PtrMesh mesh = solverTwo->_meshContexts[0]->mesh;
    utils::Vector3D integral(0.0);
    couplingInterface.integrateVectorData(dataID, raw(integral));
    validateWithMessage(tarch::la::equals(integral, utils::Vector3D(16.0)), integral);

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
    int forcesID = couplingInterface.getDataID("Forces");
    int velocitiesID = couplingInterface.getDataID("Velocities");

    int meshID = couplingInterface.getMeshID("Square");
    double pos[2];
    // Set write positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setWritePosition(meshID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setWritePosition(meshID, pos);
    // Set read positions
    pos[0] = 0.0; pos[1] = 0.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 0.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 1.0; pos[1] = 1.0;
    couplingInterface.setReadPosition(meshID, pos);
    pos[0] = 0.0; pos[1] = 1.0;
    couplingInterface.setReadPosition(meshID, pos);

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
      mesh::PtrData localForces = dataContext->localData;
      //validateEquals(forcesID, localForces->getID());
      utils::DynVector& forceValues = localForces->values();
      assign(forceValues) = 1.0;
      iterationCount++;
      maxDt = couplingInterface.advance(maxDt);

      couplingInterface.mapReadData(meshID);
      dataContext = solverOne->_dataContexts[velocitiesID];
      validate(dataContext != NULL);
      mesh::PtrData localVelocities = dataContext->localData;
      //validateEquals(velocitiesID, localVelocities->getID());
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
    int forcesID = couplingInterface.getDataID("Forces");
    int velocitiesID = couplingInterface.getDataID("Velocities");
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
      mesh::PtrData velocities = dataContext->data;
      validateEquals(velocitiesID, velocities->getID());
      utils::DynVector& values = velocities->values();
      assign(values) = 2.0;

      maxDt = couplingInterface.advance(maxDt);

      dataContext = solverTwo->_dataContexts[forcesID];
      validate(dataContext != NULL);
      mesh::PtrData forces = dataContext->localData;
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
    int forcesID = couplingInterface.getDataID("Forces");
    int velocitiesID = couplingInterface.getDataID("Velocities");
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());

    impl::PtrParticipant solverOne = couplingInterface._impl->_participants[0];
    validateEquals(solverOne->getName(), "SolverOne");
    impl::DataContext* dataContext = solverOne->_dataContexts[velocitiesID];
    validate(dataContext != NULL);
    mesh::PtrData localVelocities = dataContext->localData;
    //validateEquals(velocitiesID, localVelocities->getID());
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
    int forcesID = couplingInterface.getDataID("Forces");
    int velocitiesID = couplingInterface.getDataID("Velocities");
    double maxDt = couplingInterface.initialize();
    validate(couplingInterface.isActionRequired(actionReadSimulationCheckpoint()));
    couplingInterface.fulfilledAction(actionReadSimulationCheckpoint());

    impl::PtrParticipant solverTwo = couplingInterface._impl->_participants[1];
    validateEquals(solverTwo->getName(), "SolverTwo");
    impl::DataContext* dataContext = solverTwo->_dataContexts[forcesID];
    validate(dataContext != NULL);
    mesh::PtrData forces = dataContext->data;
    //validateEquals(velocitiesID, localVelocities->getID());
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
  std::string meshForces = "MeshForces";
  std::string meshDispl = "MeshDisplacements";
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
    int meshForcesID = interface.getMeshID(meshForces);
    int meshDisplID = interface.getMeshID(meshDispl);
    int dataForcesID = interface.getDataID(dataForces);
    int dataDisplID = interface.getDataID(dataDispl);

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
      // Set solver mesh positions for reading and writing data with mappings
      for (size_t i=0; i < size; i++){
        interface.setWritePosition(meshForcesID, raw(positions[i] + 0.1));
        interface.setReadPosition(meshDisplID, raw(positions[i] + 0.6));
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
      interface.mapReadData(meshDisplID);
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
      interface.mapReadData(meshDisplID);
      //precicePrint("2: mapped data: " << interface._impl->_accessor->dataContext(dataDisplID).data->values());
      for (size_t i=0; i < size; i++){
        interface.readVectorData(dataDisplID, i, raw(displ));
        validateNumericalEquals(displ[0], 2.0*(positions[i][0] + 0.1));
      }
      interface.finalize();
    }
    else {
      assertion1(rank == 1, rank);
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
    int meshID = precice.getMeshID("WetSurface");
    int forcesID = precice.getDataID(precice::constants::dataForces());
    int displacementsID = precice.getDataID(precice::constants::dataDisplacements());
    int oldDisplacementsID = precice.getDataID("OldDisplacements");
    validateEquals(precice.getDimensions(), 3);
    foreach(Vector3D& coord, coords){
      precice.setReadPosition(meshID, raw(coord));
      precice.setWritePosition(meshID, raw(coord));
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
      precice.mapReadData(meshID);
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
    int meshID = precice.getMeshID("WetSurface");
    foreach (Vector3D& coord, coords){
      precice.setMeshVertex(meshID, raw(coord));
    }
    for(int i=0; i < slices-1; i++){
      // Build cylinder/channel geometry
      precice.setMeshTriangleWithEdges(meshID, i*4, (i*4)+1, (i+1)*4);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4, ((i+1)*4)+1, (i*4)+1);
      precice.setMeshTriangleWithEdges(meshID, i*4+1, (i*4)+2, (i+1)*4+1);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+1, ((i+1)*4)+2, (i*4)+2);
      precice.setMeshTriangleWithEdges(meshID, i*4+2, (i*4)+3, (i+1)*4+2);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+2, ((i+1)*4)+3, (i*4)+3);
      precice.setMeshTriangleWithEdges(meshID, i*4+3, (i*4), (i+1)*4+3);
      precice.setMeshTriangleWithEdges(meshID, (i+1)*4+3, (i+1)*4, i*4);
    }
    double dt = precice.initialize();
    while(precice.isCouplingOngoing()){
      precice.advance(dt);
    }
    precice.finalize();
  }
}

#endif // defined( not PRECICE_NO_MPI )

}} // namespace precice, tests

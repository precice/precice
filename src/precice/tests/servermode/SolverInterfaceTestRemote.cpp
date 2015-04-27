// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverInterfaceTestRemote.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/RequestManager.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/Globals.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "query/tests/GeometryTestScenarios.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerIntegrationTest(precice::tests::SolverInterfaceTestRemote)

namespace precice {
namespace tests {

tarch::logging::Log SolverInterfaceTestRemote::
   _log("precice::tests::SolverInterfaceTestRemote");

SolverInterfaceTestRemote:: SolverInterfaceTestRemote()
:
  tarch::tests::TestCase("tests::SolverInterfaceTestRemote"),
  _pathToTests()
{}

SolverInterfaceTestRemote:: ~SolverInterfaceTestRemote()
{}

void SolverInterfaceTestRemote:: setUp()
{
  _pathToTests = utils::Globals::getPathToSources() + "/precice/tests/servermode/";
}

void SolverInterfaceTestRemote:: run()
{
# ifndef PRECICE_NO_MPI
  preciceTrace("run()");
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() >= 2){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1;
    Par::Communicator comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 1){
      Par::setGlobalCommunicator(comm);
      testMethod(testGeometryMode);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  if (Par::getCommunicatorSize() >= 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2;
    Par::Communicator comm = Par::getRestrictedCommunicator(ranksWanted);
    if ( Par::getProcessRank() <= 2 ){
      Par::setGlobalCommunicator(comm);
      testMethod(testCouplingModeWithOneServer);
      testMethod(testGeometryModeParallel);
      testMethod(testGeometryModeParallelStationaryMapping);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  if (Par::getCommunicatorSize() >= 4){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2, 3;
    Par::Communicator comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod(testCouplingModeParallelWithOneServer);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

void SolverInterfaceTestRemote:: configureSolverInterface
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

void SolverInterfaceTestRemote:: testGeometryMode()
{
  preciceTrace("testGeometryMode()");
  using namespace tarch::la;
  using utils::DynVector;
  for (int dim=2; dim <= 3; dim++){
    std::string configFilename;
    if (dim == 2){
      configFilename = _pathToTests + "geomode-2D.xml";
    }
    else {
      configFilename = _pathToTests + "geomode-3D.xml";
    }
    int rank = utils::Parallel::getProcessRank();
    if (rank == 0){
      SolverInterface interface("TestAccessor", rank, 2);
      configureSolverInterface(configFilename, interface);
      validateEquals(interface.getDimensions(), dim);

      int meshIDScalar = interface.getMeshID ( "AccessorMeshScalar" );
      int meshIDVector = interface.getMeshID ( "AccessorMeshVector" );
      DynVector pos(dim, 4.0);
      int posIndex1 = interface.setMeshVertex ( meshIDScalar, raw(pos) );
      assign(pos) = 3.0;
      int posIndex2 = interface.setMeshVertex ( meshIDScalar, raw(pos) );
      assign(pos) = 3.0;
      int posIndex3 = interface.setMeshVertex ( meshIDVector, raw(pos) );

      interface.initialize();
      interface.initializeData(); // is skipped due to geometry mode

      std::set<int> ids = interface.getMeshIDs();
      typedef query::tests::GeometryTestScenarios GeoTests;
      GeoTests geoTests;

      // Test inquireClosestMesh()
      const GeoTests::PointQueryScenario& pointScen = geoTests.pointQueryScenario(dim);
      std::list<DynVector>::const_iterator coordIter = pointScen.queryCoords.begin();
      std::list<double>::const_iterator distIter = pointScen.validDistances.begin();
      std::list<DynVector>::const_iterator distVectorIter = pointScen.validDistanceVectors.begin();
      DynVector distanceVec(dim);
      while (coordIter != pointScen.queryCoords.end()){
        ClosestMesh closest = interface.inquireClosestMesh(raw(*coordIter), ids);
        for(int i=0; i<dim; i++) distanceVec[i] = closest.distanceVector()[i];
        //validate(equals(*distVectorIter, distanceVec));
        //validate(equals(*distIter, closest.distance()));
        coordIter++;
        distIter++;
        distVectorIter++;
      }

      // Test inquirePosition()
      const GeoTests::PositionQueryScenario& posScen = geoTests.positionQueryScenario(dim);
      coordIter = posScen.queryCoords.begin();
      std::list<int>::const_iterator posIter = posScen.validPositions.begin();
      while ( coordIter != posScen.queryCoords.end() ){
        int position = interface.inquirePosition ( raw(*coordIter), ids );
        validateEquals ( position, *posIter );
        coordIter ++;
        posIter ++;
      }

      // Test inquireVoxelPosition()
      const GeoTests::VoxelQueryScenario& voxelScen = geoTests.voxelQueryScenario(dim);
      std::list<DynVector>::const_iterator centerIter = voxelScen.queryCenters.begin();
      std::list<DynVector>::const_iterator hIter = voxelScen.queryHalflengths.begin();
      std::list<bool>::const_iterator includeBoundsIter = voxelScen.includeBoundaries.begin();
      posIter = voxelScen.validPositions.begin();
      while ( centerIter != voxelScen.queryCenters.end() ){
        VoxelPosition pos = interface.inquireVoxelPosition ( raw(*centerIter),
                            raw(*hIter), *includeBoundsIter, ids );
        validateEquals ( pos.position(), *posIter );
        centerIter ++;
        hIter ++;
        includeBoundsIter ++;
        posIter ++;
      }



      // Test write data

      int dataID = interface.getDataID ( "ScalarData", meshIDScalar );
      double value = 1.0;
      interface.writeScalarData ( dataID, posIndex1, value );
      value = 2.0;
      interface.writeScalarData ( dataID, posIndex2, value );

      // Test read data (not really good test...)
      dataID = interface.getDataID ( "VectorData", meshIDVector );
      DynVector readValue(dim, 2.0);
      interface.readVectorData ( dataID, posIndex3, raw(readValue) );
      validate ( equals(readValue, DynVector(dim,0.0)) );

      // Test exporting mesh
      interface.exportMesh ( "remote" );

      interface.advance(1.0);

      interface.finalize();
    }
    else {
      assertion1 ( rank == 1, rank );
      bool isServer = true;
      impl::SolverInterfaceImpl server ( "TestAccessor", rank, 2, isServer );
      // Perform manual configuration without overwritting logging config
      mesh::Mesh::resetGeometryIDsGlobally();
      mesh::Data::resetDataCount();
      impl::Participant::resetParticipantCount();
      config::Configuration config;
      utils::configure(config.getXMLTag(), configFilename);
      server.configure(config.getSolverInterfaceConfiguration());

      validateEquals(server.getDimensions(), dim);
      server.runServer();
    }
  }
}

void SolverInterfaceTestRemote:: testGeometryModeParallel()
{
  preciceTrace("testGeometryModeParalell()");
  using namespace tarch::la;
  using utils::DynVector;
  for ( int dim=2; dim <= 3; dim++ ){
    std::string configFilename;
    if (dim == 2){
      configFilename = _pathToTests + "geomode-2D.xml";
    }
    else {
      configFilename = _pathToTests + "geomode-3D.xml";
    }
    int rank = utils::Parallel::getProcessRank();
    if ( (rank == 0) || (rank == 1) ){
      SolverInterface interface ( "TestAccessor", rank, 2 );
      configureSolverInterface ( configFilename, interface );
      validateEquals ( interface.getDimensions(), dim );
      interface.initialize();
      interface.initializeData(); // is skipped due to geometry mode

      int meshID = interface.getMeshID("CuboidMesh");
      std::set<int> ids;
      ids.insert(meshID);
      typedef query::tests::GeometryTestScenarios GeoTests;
      GeoTests geoTests;

      // Test inquireClosestMesh()
      const GeoTests::PointQueryScenario& pointScen = geoTests.pointQueryScenario(dim);
      std::list<DynVector>::const_iterator coordIter = pointScen.queryCoords.begin();
      std::list<double>::const_iterator distIter = pointScen.validDistances.begin();
      std::list<DynVector>::const_iterator distVectorIter = pointScen.validDistanceVectors.begin();
      DynVector distanceVec(dim);
      while ( coordIter != pointScen.queryCoords.end() ){
        ClosestMesh closest = interface.inquireClosestMesh ( raw(*coordIter), ids );
        for(int i=0; i<dim; i++) distanceVec[i] = closest.distanceVector()[i];
        validate ( equals(*distVectorIter, distanceVec) );
        validate ( equals(*distIter, closest.distance()) );
        coordIter ++;
        distIter ++;
        distVectorIter ++;
      }

      // Test inquirePosition()
      const GeoTests::PositionQueryScenario & posScen = geoTests.positionQueryScenario(dim);
      coordIter = posScen.queryCoords.begin();
      std::list<int>::const_iterator posIter = posScen.validPositions.begin();
      while ( coordIter != posScen.queryCoords.end() ){
        int position = interface.inquirePosition ( raw(*coordIter), ids );
        validateEquals ( position, *posIter );
        coordIter ++;
        posIter ++;
      }

      // Test inquireVoxelPosition()
      const GeoTests::VoxelQueryScenario& voxelScen = geoTests.voxelQueryScenario(dim);
      std::list<DynVector>::const_iterator centerIter = voxelScen.queryCenters.begin();
      std::list<DynVector>::const_iterator hIter = voxelScen.queryHalflengths.begin();
      std::list<bool>::const_iterator includeBoundsIter = voxelScen.includeBoundaries.begin();
      posIter = voxelScen.validPositions.begin();
      while ( centerIter != voxelScen.queryCenters.end() ){
        VoxelPosition pos = interface.inquireVoxelPosition ( raw(*centerIter),
                            raw(*hIter), *includeBoundsIter, ids );
        validateEquals ( pos.position(), *posIter );
        centerIter ++;
        hIter ++;
        includeBoundsIter ++;
        posIter ++;
      }

      int meshIDScalar = interface.getMeshID ( "AccessorMeshScalar" );
      int meshIDVector = interface.getMeshID ( "AccessorMeshVector" );

      // Test write data
      DynVector pos(dim, 0.0);
      int posIndex = interface.setMeshVertex(meshIDScalar, raw(pos));
      int dataID = interface.getDataID("ScalarData", meshIDScalar);
      double value = 1.0;
      interface.writeScalarData(dataID, posIndex, value);

      assign(pos) = 1.0;
      posIndex = interface.setMeshVertex(meshIDScalar, raw(pos));
      value = 2.0;
      interface.writeScalarData(dataID, posIndex, value);

      // Let the written data of both processes be accumulated
      interface._impl->_requestManager->requestPing();
      utils::Parallel::synchronizeLocalProcesses();

      // Test read data (not really good test...)
      assign(pos) = 0.0;
      posIndex = interface.setMeshVertex(meshIDVector, raw(pos));
      dataID = interface.getDataID("VectorData", meshIDVector);
      DynVector readValue(dim, 4.0);
      interface.readVectorData(dataID, posIndex, raw(readValue));
      validate(equals(readValue, DynVector(dim,0.0)));

      // Test exporting mesh
      std::ostringstream filename;
      filename << "remote-parallel-rank" << rank;
      interface.exportMesh(filename.str());

      interface.advance(1.0);

      interface.finalize();
    }
    else {
      assertion1(rank == 2, rank);
      bool isServer = true;
      impl::SolverInterfaceImpl server("TestAccessor", 0, 1, isServer);
      // Perform manual configuration without overwritting logging config
      mesh::Mesh::resetGeometryIDsGlobally();
      mesh::Data::resetDataCount();
      impl::Participant::resetParticipantCount();
      config::Configuration config;
      utils::configure(config.getXMLTag(), configFilename);
      server.configure(config.getSolverInterfaceConfiguration());

      validateEquals(server.getDimensions(), dim);
      server.runServer();
    }
  }
}

void SolverInterfaceTestRemote:: testGeometryModeParallelStationaryMapping()
{
  preciceTrace("testGeometryModeParallelStationaryMapping()");
  using namespace tarch::la;
  using utils::DynVector;
  for ( int dim=2; dim <= 2; dim++ ){ // LIMITED TO 2D!!!!!!
    std::string configFilename;
    if (dim == 2){
      configFilename = _pathToTests + "geomode-stationary-mapping-2D.xml";
    }
    else {
      configFilename = _pathToTests + "geomode-stationary-mapping-3D.xml";
    }
    int rank = utils::Parallel::getProcessRank();
    if ((rank == 0) || (rank == 1)){
      SolverInterface interface("Accessor", rank, 2);
      configureSolverInterface(configFilename, interface);
      validateEquals(interface.getDimensions(), dim);
      interface.initialize();
      interface.initializeData(); // is skipped due to geometry mode

      // Test write data
      DynVector pos(dim, 0.0);
      int meshIDVector = interface.getMeshID ( "AccessorMeshVector" );
      int indices[4];

      //int scalarDataID = interface.getDataID("ScalarData");
      int vectorDataID = interface.getDataID("VectorData", meshIDVector);
      //double scalarValues[] = {1.0, 3.0, 4.0, 2.0};

      if (rank == 0){
        pos[0] = 0.0; pos[1] = 0.0;
        indices[0] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 1.0; pos[1] = 0.0;
        indices[2] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 1.0; pos[1] = 1.0;
        indices[3] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 0.0; pos[1] = 1.0;
        indices[1] = interface.setMeshVertex(meshIDVector, raw(pos));
      }
      else {
        pos[0] = 0.5; pos[1] = 0.0;
        indices[0] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 1.0; pos[1] = 0.5;
        indices[2] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 0.5; pos[1] = 1.0;
        indices[3] = interface.setMeshVertex(meshIDVector, raw(pos));
        pos[0] = 0.0; pos[1] = 0.5;
        indices[1] = interface.setMeshVertex(meshIDVector, raw(pos));
      }

      double vectorValues[] = {1.0, 1.0, 3.0, 3.0, 4.0, 4.0, 2.0, 2.0};

      // IMPLEMENT WRITE BLOCK SCALAR DATA
      interface.writeBlockVectorData(vectorDataID, 4, indices, vectorValues);
      interface.mapWriteDataFrom(meshIDVector);

      // Let the written data of both processes be accumulated
//      interface._impl->_accessor->getClientServerCommunication()->send(
//          impl::SolverInterfaceImpl::REQUEST_PING, 0);
//      bool ping = false;
//      interface._impl->_accessor->getClientServerCommunication()->receive(ping, 0);
//      utils::Parallel::synchronizeLocalProcesses();
      interface.exportMesh("testGeometryModeParallelStationaryMapping");

      interface.finalize();
    }
    else {
      assertion1(rank == 2, rank);
      bool isServer = true;
      impl::SolverInterfaceImpl server("Accessor", 0, 1, isServer);
      // Perform manual configuration without overwritting logging config
      mesh::Mesh::resetGeometryIDsGlobally();
      mesh::Data::resetDataCount();
      impl::Participant::resetParticipantCount();
      config::Configuration config;
      utils::configure(config.getXMLTag(), configFilename);
      server.configure(config.getSolverInterfaceConfiguration());
      validateEquals(server.getDimensions(), dim);
      server.runServer();
    }
  }
}

void SolverInterfaceTestRemote:: testCouplingModeWithOneServer()
{
  preciceTrace( "testCouplingModeWithOneServer()" );
  int rank = utils::Parallel::getProcessRank();
  std::string configFile = _pathToTests + "cplmode-1.xml";
  if ( rank == 0 ){
    SolverInterface interface("ParticipantA", 0, 1);
    configureSolverInterface ( configFile, interface );
    double time = 0.0;
    int timesteps = 0;
    double dt = interface.initialize();
    while ( interface.isCouplingOngoing() ){
      time += dt;
      timesteps++;
      dt = interface.advance(dt);
    }
    interface.finalize();
    validateNumericalEquals ( time, 5.0 );
    validateEquals ( timesteps, 5 );
  }
  else if ( rank == 1 ){
    SolverInterface interface("ParticipantB", 0, 1);
    configureSolverInterface ( configFile, interface );
    double time = 0.0;
    int timesteps = 0;
    double dt = interface.initialize();
    while ( interface.isCouplingOngoing() ){
      time += dt;
      timesteps++;
      dt = interface.advance(dt);
    }
    interface.finalize();
    validateNumericalEquals ( time, 5.0 );
    validateEquals ( timesteps, 5 );
  }
  else {
    assertion1 (rank == 2, rank);
    bool isServer = true;
    impl::SolverInterfaceImpl server("ParticipantB", 0, 1, isServer);

    // Perform manual configuration without overwritting logging config
    mesh::Mesh::resetGeometryIDsGlobally();
    mesh::Data::resetDataCount();
    impl::Participant::resetParticipantCount();
    config::Configuration config;
    utils::configure ( config.getXMLTag(), configFile );
    //validate ( config.isValid() );
    server.configure ( config.getSolverInterfaceConfiguration() );

    server.runServer();
  }
}

void SolverInterfaceTestRemote:: testCouplingModeParallelWithOneServer()
{
  preciceTrace( "testCouplingModeParallelWithOneServer()" );
  using namespace tarch::la;
  int rank = utils::Parallel::getProcessRank();
  std::string configFile = _pathToTests + "cplmode-1.xml";
  if (rank == 0){
    SolverInterface interface("ParticipantA", 0, 1);
    configureSolverInterface(configFile, interface);
    double time = 0.0;
    int timesteps = 0;
    double dt = interface.initialize();
    MeshHandle handle = interface.getMeshHandle("Mesh");
    VertexHandle vertices = handle.vertices();
    int meshID = interface.getMeshID("Mesh");
    int scalarDataID = interface.getDataID("ScalarData", meshID);
    int vectorDataID = interface.getDataID("VectorData", meshID);
    int dataSize = 4;
    int indices[] = {0, 1, 2, 3};
    double vectorValues[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    utils::DynVector expect(8);
    assignList(expect) = 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0;
    while (interface.isCouplingOngoing()){
      time += dt;
      timesteps++;
      foriter (VertexIterator, vertex, vertices){
        interface.writeScalarData(scalarDataID, vertex.vertexID(), 1.0);
      }
      dt = interface.advance(dt);
      interface.readBlockVectorData(vectorDataID, dataSize, indices, vectorValues);
      validateWithParams2(equals(wrap<8>(vectorValues), expect),
                          wrap<8>(vectorValues), expect);
      assign(wrap<8>(vectorValues)) = 0.0;
    }
    interface.finalize();
    validateNumericalEquals(time, 5.0);
    validateEquals(timesteps, 5);
  }
  else if ((rank == 1) || (rank == 2)){
    SolverInterface interface("ParticipantB", rank-1, 2);
    configureSolverInterface(configFile, interface);
    double time = 0.0;
    int timesteps = 0;
    double dt = interface.initialize();
    int meshID = interface.getMeshID("Mesh");
    int scalarDataID = interface.getDataID("ScalarData", meshID);
    int vectorDataID = interface.getDataID("VectorData", meshID);
    int dataSize = 4;
    int indices[] = {0, 1, 2, 3};
    double vectorValues[] = {1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
    while (interface.isCouplingOngoing()){
      double value = 0.0;
      interface.readScalarData(scalarDataID, 0, value);
      validateWithParams1(tarch::la::equals(value, 1.0), value);
      interface.writeBlockVectorData(vectorDataID, dataSize, indices, vectorValues);
      time += dt;
      timesteps++;
      dt = interface.advance(dt);
    }
    interface.finalize();
    validateNumericalEquals(time, 5.0);
    validateEquals(timesteps, 5);
  }
  else {
    assertion1(rank == 3, rank);
    bool isServer = true;
    impl::SolverInterfaceImpl server("ParticipantB", 0, 1, isServer);

    // Perform manual configuration without overwritting logging config
    mesh::Mesh::resetGeometryIDsGlobally();
    mesh::Data::resetDataCount();
    impl::Participant::resetParticipantCount();
    config::Configuration config;
    utils::configure(config.getXMLTag(), configFile);
    //validate ( config.isValid() );
    server.configure(config.getSolverInterfaceConfiguration());
    server.runServer();
  }
}

}} // namespace precice, tests

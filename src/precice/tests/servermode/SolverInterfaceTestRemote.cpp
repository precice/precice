#include "SolverInterfaceTestRemote.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/RequestManager.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/Globals.hpp"
#include "com/MPIDirectCommunication.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerIntegrationTest(precice::tests::SolverInterfaceTestRemote)

namespace precice {
namespace tests {

logging::Logger SolverInterfaceTestRemote::
   _log("tests::SolverInterfaceTestRemote");

SolverInterfaceTestRemote:: SolverInterfaceTestRemote()
:
  tarch::tests::TestCase("tests::SolverInterfaceTestRemote"),
  _pathToTests()
{}

SolverInterfaceTestRemote:: ~SolverInterfaceTestRemote()
{}

void SolverInterfaceTestRemote:: setUp()
{
  _pathToTests = utils::getPathToSources() + "/precice/tests/servermode/";
}

void SolverInterfaceTestRemote:: run()
{
# ifndef PRECICE_NO_MPI
  TRACE();
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() >= 3){
    Par::Communicator comm = Par::getRestrictedCommunicator({0, 1, 2});
    if ( Par::getProcessRank() <= 2 ){
      Par::setGlobalCommunicator(comm);
      testMethod(testCouplingModeWithOneServer);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
  if (Par::getCommunicatorSize() >= 4){
    Par::Communicator comm = Par::getRestrictedCommunicator({0, 1, 2, 3});
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
  TRACE(configFilename);
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  config::Configuration config;
  xml::configure(config.getXMLTag(), configFilename);
  //validate ( config.isValid() );
  interface._impl->configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceTestRemote:: testCouplingModeWithOneServer()
{
  TRACE();
  int rank = utils::Parallel::getProcessRank();
  std::string configFile = _pathToTests + "cplmode-1.xml";
  if ( rank == 0 ){
    SolverInterface interface("ParticipantA", 0, 1);
    configureSolverInterface ( configFile, interface );
    double time = 0.0;
    int timesteps = 0;

    int meshID = interface.getMeshID("Mesh");
    interface.setMeshVertex(meshID, Eigen::Vector2d(0.0,0.0).data());
    interface.setMeshVertex(meshID, Eigen::Vector2d(1.0,0.0).data());

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
    assertion (rank == 2, rank);
    bool isServer = true;
    impl::SolverInterfaceImpl server("ParticipantB", 0, 1, isServer);

    // Perform manual configuration without overwritting logging config
    mesh::Mesh::resetGeometryIDsGlobally();
    mesh::Data::resetDataCount();
    impl::Participant::resetParticipantCount();
    config::Configuration config;
    xml::configure ( config.getXMLTag(), configFile );
    //validate ( config.isValid() );
    server.configure ( config.getSolverInterfaceConfiguration() );

    server.runServer();
  }
}

void SolverInterfaceTestRemote:: testCouplingModeParallelWithOneServer()
{
  TRACE();
  int rank = utils::Parallel::getProcessRank();
  std::string configFile = _pathToTests + "cplmode-1.xml";
  if (rank == 0){
    SolverInterface interface("ParticipantA", 0, 1);
    configureSolverInterface(configFile, interface);
    double time = 0.0;
    int timesteps = 0;

    int meshID = interface.getMeshID("Mesh");
    int scalarDataID = interface.getDataID("ScalarData", meshID);
    int vectorDataID = interface.getDataID("VectorData", meshID);
    interface.setMeshVertex(meshID, Eigen::Vector2d(1.0,0.0).data());
    interface.setMeshVertex(meshID, Eigen::Vector2d(0.0,-1.0).data());
    interface.setMeshVertex(meshID, Eigen::Vector2d(-1.0,0.0).data());
    interface.setMeshVertex(meshID, Eigen::Vector2d(0.0,1.0).data());

    double dt = interface.initialize();
    MeshHandle handle = interface.getMeshHandle("Mesh");
    VertexHandle vertices = handle.vertices();
    int dataSize = 4;
    int indices[] = {0, 1, 2, 3};
    double vectorValues[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Eigen::Matrix<double, 8, 1> expect;
    expect << 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0;
    while (interface.isCouplingOngoing()){
      time += dt;
      timesteps++;
      for (VertexIterator vertex : vertices) {
        interface.writeScalarData(scalarDataID, vertex.vertexID(), 1.0);
      }
      dt = interface.advance(dt);
      interface.readBlockVectorData(vectorDataID, dataSize, indices, vectorValues);
      validate(math::equals(Eigen::Map<Eigen::Matrix<double, 8, 1>>(vectorValues), expect));              
      Eigen::Map<Eigen::Matrix<double, 8, 1>>(vectorValues).setConstant(0);
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
      validateWithParams1(math::equals(value, 1.0), value);
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
    assertion(rank == 3, rank);
    bool isServer = true;
    impl::SolverInterfaceImpl server("ParticipantB", 0, 1, isServer);

    // Perform manual configuration without overwritting logging config
    mesh::Mesh::resetGeometryIDsGlobally();
    mesh::Data::resetDataCount();
    impl::Participant::resetParticipantCount();
    config::Configuration config;
    xml::configure(config.getXMLTag(), configFile);
    //validate ( config.isValid() );
    server.configure(config.getSolverInterfaceConfiguration());
    server.runServer();
  }
}

}} // namespace precice, tests

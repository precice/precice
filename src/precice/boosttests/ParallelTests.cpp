#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"

#include "utils/Parallel.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include <fstream>

using namespace precice;

namespace precice {
extern bool testMode;
}


void reset ()
{
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  precice::utils::EventRegistry::clear();
  utils::MasterSlave::reset();
}

struct ParallelTestFixture {
  ParallelTestFixture()
  {
    utils::Parallel::restrictGlobalCommunicator({0,1,2,3});
    assertion(utils::Parallel::getCommunicatorSize() == 4);
    precice::testMode = true;
  }

  ~ParallelTestFixture()
  {
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    reset();
  }
};



BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Parallel, ParallelTestFixture)


BOOST_AUTO_TEST_CASE(TestMasterSlaveSetup, * testing::OnSize(4))
{
  SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 4 );
  std::string configFilename = utils::getPathToSources() + "/precice/boosttests/config1.xml";
  config::Configuration config;
  utils::configure(config.getXMLTag(), configFilename);
  interface._impl->configure(config.getSolverInterfaceConfiguration());

  BOOST_TEST ( interface.getDimensions() == 3 );

  if(utils::Parallel::getProcessRank()==0){
    BOOST_TEST(utils::MasterSlave::_masterMode == true);
    BOOST_TEST(utils::MasterSlave::_slaveMode == false);
  }
  else {
    BOOST_TEST(utils::MasterSlave::_masterMode == false);
    BOOST_TEST(utils::MasterSlave::_slaveMode == true);
  }

  BOOST_TEST(utils::MasterSlave::_rank == utils::Parallel::getProcessRank());
  BOOST_TEST(utils::MasterSlave::_size == 4);
  BOOST_TEST(utils::MasterSlave::_communication.use_count()>0);
  BOOST_TEST(utils::MasterSlave::_communication->isConnected());


  //necessary as this test does not call finalize
  utils::MasterSlave::_communication = nullptr;
  utils::Parallel::clearGroups();
}

BOOST_AUTO_TEST_CASE(TestFinalize, * testing::OnSize(4))
{
  std::string configFilename = utils::getPathToSources() + "/precice/boosttests/config1.xml";
  config::Configuration config;
  utils::configure(config.getXMLTag(), configFilename);
  if(utils::Parallel::getProcessRank()<=1){
    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 2 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshOne");
    double xCoord = 0.0 + utils::Parallel::getProcessRank();
    interface.setMeshVertex(meshID, Eigen::Vector3d(xCoord,0.0,0.0).data());
    interface.initialize();
    BOOST_TEST(interface.getMeshHandle("MeshOne").vertices().size()==1);
    BOOST_TEST(interface.getMeshHandle("MeshTwo").vertices().size()==1);
    interface.finalize();
  }
  else {
    SolverInterface interface ( "SolverTwo", utils::Parallel::getProcessRank()-2, 2 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshTwo");
    double xCoord = -2.0 + utils::Parallel::getProcessRank();
    interface.setMeshVertex(meshID, Eigen::Vector3d(xCoord,0.0,0.0).data());
    interface.initialize();
    BOOST_TEST(interface.getMeshHandle("MeshTwo").vertices().size()==1);
    interface.finalize();
  }
}

#ifndef PRECICE_NO_PETSC
BOOST_AUTO_TEST_CASE(GlobalRBFPartitioning, * testing::OnSize(4))
{
  std::string configFilename = utils::getPathToSources() + "/precice/boosttests/globalRBFPartitioning.xml";
  config::Configuration config;

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    utils::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshOne");
    int dataID = interface.getDataID("Data2", meshID);

    int vertexIDs[2];
    double xCoord = utils::Parallel::getProcessRank() * 0.4;
    double positions[4] = {xCoord,0.0,xCoord+0.2,0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    double values[2];
    interface.advance(1.0);
    interface.readBlockScalarData(dataID, 2, vertexIDs, values);
//    std::cout << utils::Parallel::getProcessRank() <<": " << values << std::endl;
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "SolverTwo" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    utils::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverTwo", 0, 1 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshTwo");
    int vertexIDs[6];
    double positions[12] = {0.0,0.0,0.2,0.0,0.4,0.0,0.6,0.0,0.8,0.0,1.0,0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    int dataID = interface.getDataID("Data2", meshID);
    double values[6] = {1.0,2.0,3.0,4.0,5.0,6.0};
    interface.writeBlockScalarData(dataID, 6, vertexIDs, values);
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(LocalRBFPartitioning, * testing::OnSize(4))
{
  std::string configFilename = utils::getPathToSources() + "/precice/boosttests/localRBFPartitioning.xml";
  config::Configuration config;

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    utils::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshOne");
    int dataID = interface.getDataID("Data2", meshID);

    int vertexIDs[2];
    double xCoord = utils::Parallel::getProcessRank() * 0.4;
    double positions[4] = {xCoord,0.0,xCoord+0.2,0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    double values[2];
    interface.advance(1.0);
    interface.readBlockScalarData(dataID, 2, vertexIDs, values);
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "SolverTwo" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    utils::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverTwo", 0, 1 );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("MeshTwo");
    int vertexIDs[6];
    double positions[12] = {0.0,0.0,0.2,0.0,0.4,0.0,0.6,0.0,0.8,0.0,1.0,0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    int dataID = interface.getDataID("Data2", meshID);
    double values[6] = {1.0,2.0,3.0,4.0,5.0,6.0};
    interface.writeBlockScalarData(dataID, 6, vertexIDs, values);
    interface.advance(1.0);
    interface.finalize();
  }
}

#endif // PRECICE_NO_PETSC

/// tests for various QN settings if correct number of iterations is returned
BOOST_AUTO_TEST_CASE(TestQN, * testing::OnSize(4))
{
  int numberOfTests = 3;
  std::vector<std::string> configs;
  configs.resize(numberOfTests);
  configs[0] = utils::getPathToSources() + "/precice/boosttests/QN1.xml";
  configs[1] = utils::getPathToSources() + "/precice/boosttests/QN2.xml";
  configs[2] = utils::getPathToSources() + "/precice/boosttests/QN3.xml";

  int correctIterations[3] = {29, 17, 15};

  std::string solverName, meshName, writeDataName, readDataName;
  int rank, size;

  if(utils::Parallel::getProcessRank()==0){
    solverName = "SolverOne";
    meshName = "MeshOne";
    writeDataName = "Data1";
    readDataName = "Data2";
    rank = 0;
    size = 1;
  }
  else {
    solverName = "SolverTwo";
    meshName = "MeshTwo";
    writeDataName = "Data2";
    readDataName = "Data1";
    rank = utils::Parallel::getProcessRank() - 1;
    size = 3;
  }

  utils::Parallel::splitCommunicator( solverName );
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
  assertion(utils::Parallel::getCommunicatorSize() == size);
  utils::Parallel::clearGroups();

  for(int k=0; k<numberOfTests; k++){
    reset();
    config::Configuration config;
    std::string configFilename = configs[k];

    utils::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( solverName, rank, size );
    interface._impl->configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID(meshName);
    int writeDataID = interface.getDataID(writeDataName, meshID);
    int readDataID = interface.getDataID(readDataName, meshID);

    int vertexIDs[4];

    if(utils::Parallel::getProcessRank()==0){
      double positions[8] = {2.0, 0.0, 2.0, 0.5, 2.0, 1.0, 2.5, 1.0};
      interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    }
    else if(utils::Parallel::getProcessRank()==1){
      double positions[8] = {2.0, 0.1, 2.0, 0.25, 2.0, 0.4, 2.0, 0.5};
      interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    }
    else if(utils::Parallel::getProcessRank()==2){
      double positions[8] = {2.0, 0.6, 2.0, 0.75, 2.0, 0.9, 2.0, 1.0};
      interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    }
    else if(utils::Parallel::getProcessRank()==3){
      double positions[8] = {2.1, 1.0, 2.25, 1.0, 2.4, 1.0, 2.5, 1.0};
      interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    }

    interface.initialize();
    double inValues[4] = {0.0, 0.0, 0.0, 0.0};
    double outValues[4] = {0.0, 0.0, 0.0, 0.0};

    int iterations = 0;

    while(interface.isCouplingOngoing()){
      if(interface.isActionRequired(precice::constants::actionWriteIterationCheckpoint())){
        interface.fulfilledAction(precice::constants::actionWriteIterationCheckpoint());
      }

      if(utils::Parallel::getProcessRank()==0){
        for(int i=0; i<4; i++){
          outValues[i] = inValues[i]*inValues[i] - 30.0;
        }
      }
      else {
        for(int i=0; i<4; i++){
          outValues[i] = inValues[(i+1)%4]*inValues[(i+2)%4] - 2.0;
        }
      }

      interface.writeBlockScalarData(writeDataID, 4, vertexIDs, outValues);
      interface.advance(1.0);
      interface.readBlockScalarData(readDataID, 4, vertexIDs, inValues);

      if(interface.isActionRequired(precice::constants::actionReadIterationCheckpoint())){
        interface.fulfilledAction(precice::constants::actionReadIterationCheckpoint());
        iterations++;
      }
    }
    interface.finalize();
    BOOST_TEST(iterations == correctIterations[k]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "utils/Parallel.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/Parallel.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Event.hpp"

using namespace precice;


void reset ()
{
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  utils::MasterSlave::reset();
}

struct ParallelTestFixture : testing::WhiteboxFixture {

  std::string _pathToTests;

  ParallelTestFixture()
  {
    reset();
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
    utils::Parallel::restrictGlobalCommunicator({0,1,2,3});
    assertion(utils::Parallel::getCommunicatorSize() == 4);
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
  std::string configFilename = _pathToTests + "config1.xml";
  config::Configuration config;
  xml::configure(config.getXMLTag(), configFilename);
  impl(interface).configure(config.getSolverInterfaceConfiguration());

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
  std::string configFilename = _pathToTests + "config1.xml";
  config::Configuration config;
  xml::configure(config.getXMLTag(), configFilename);
  if(utils::Parallel::getProcessRank()<=1){
    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 2 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
  std::string configFilename = _pathToTests + "globalRBFPartitioning.xml";
  config::Configuration config;
  utils::MasterSlave::_rank = utils::Parallel::getProcessRank();
  
  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator()); //needed since this test uses PETSc
    assertion(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    xml::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
//    std::cout << utils::Parallel::getProcessRank() <<": " << values << '\n';
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "SolverTwo" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    xml::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverTwo", 0, 1 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
  std::string configFilename = _pathToTests + "localRBFPartitioning.xml";
  config::Configuration config;

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    xml::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
    xml::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( "SolverTwo", 0, 1 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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

/// This testcase is based on a bug reported by Thorsten for acoustic FASTEST-Ateles coupling
BOOST_AUTO_TEST_CASE(CouplingOnLine, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "line-coupling.xml";
  config::Configuration config;

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "Ateles" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    xml::configure(config.getXMLTag(), configFilename);
    SolverInterface interface ( "Ateles", utils::Parallel::getProcessRank(), 3 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("Ateles_Mesh");

    int vertexIDs[4];
    double offset = utils::Parallel::getProcessRank() * 0.4;
    double xCoord = 0.0;
    double yCoord = 1.0;
    double positions[12] = {xCoord, yCoord, 0.1 + offset,
                            xCoord, yCoord, 0.2 + offset,
                            xCoord, yCoord, 0.3 + offset,
                            xCoord, yCoord, 0.4 + offset};
    interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "FASTEST" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    assertion(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    xml::configure(config.getXMLTag(), configFilename);
    SolverInterface interface ( "FASTEST", 0, 1 );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
    int meshID = interface.getMeshID("FASTEST_Mesh");
    int vertexIDs[10];
    double xCoord = -0.0001;
    double yCoord = 1.00001;
    double positions[30] = {xCoord, yCoord, 0.12,
                            xCoord, yCoord, 0.24,
                            xCoord, yCoord, 0.36,
                            xCoord, yCoord, 0.48,
                            xCoord, yCoord, 0.60,
                            xCoord, yCoord, 0.72,
                            xCoord, yCoord, 0.84,
                            xCoord, yCoord, 0.96,
                            xCoord, yCoord, 1.08,
                            xCoord, yCoord, 1.2};
    interface.setMeshVertices(meshID, 10, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}


/// tests for various QN settings if correct number of iterations is returned
BOOST_AUTO_TEST_CASE(TestQN, * testing::OnSize(4))
{
  int numberOfTests = 3;
  std::vector<std::string> configs;
  configs.resize(numberOfTests);
  configs[0] = _pathToTests + "QN1.xml";
  configs[1] = _pathToTests + "QN2.xml";
  configs[2] = _pathToTests + "QN3.xml";

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

    xml::configure(config.getXMLTag(), configFilename);

    SolverInterface interface ( solverName, rank, size );
    impl(interface).configure(config.getSolverInterfaceConfiguration());
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
    // Depending on the hardware, QN can be slighly faster or slower leading to an iteration more or less.
    BOOST_TEST(iterations <= correctIterations[k] + 1);
    BOOST_TEST(iterations >= correctIterations[k] - 1);
  }
}

// This test does not restrict the communicator per participant, since otherwise MPI ports do not work for Open-MPI
/// Tests various distributed communication schemes.
BOOST_AUTO_TEST_CASE(testDistributedCommunications, * testing::OnSize(4))
{
  std::vector<std::string> fileNames({
      "point-to-point-sockets.xml",
      "point-to-point-mpi.xml",
      "gather-scatter-mpi.xml"});

  for (auto fileName : fileNames) {
    reset();

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
    config::Configuration config;
    xml::configure(config.getXMLTag(), _pathToTests + fileName);
    impl(precice).configure(config.getSolverInterfaceConfiguration());
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
          BOOST_TEST(expectedData[i+i1][d] == data[i+i1][d]);
        }
      }
    }

    precice.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

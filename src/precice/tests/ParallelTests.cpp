#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "utils/Parallel.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
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

struct ParallelTestFixture : testing::SlimConfigurator {

  std::string _pathToTests;

  ParallelTestFixture()
  {
    reset();
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
    utils::Parallel::restrictGlobalCommunicator({0,1,2,3});
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);
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
  slimConfigure(interface, configFilename);
  BOOST_TEST ( interface.getDimensions() == 3 );

  if(utils::Parallel::getProcessRank()==0){
    BOOST_TEST(utils::MasterSlave::isMaster() == true);
    BOOST_TEST(utils::MasterSlave::isSlave() == false);
  }
  else {
    BOOST_TEST(utils::MasterSlave::isMaster() == false);
    BOOST_TEST(utils::MasterSlave::isSlave() == true);
  }

  BOOST_TEST(utils::MasterSlave::getRank() == utils::Parallel::getProcessRank());
  BOOST_TEST(utils::MasterSlave::getSize() == 4);
  BOOST_TEST(utils::MasterSlave::_communication.use_count()>0);
  BOOST_TEST(utils::MasterSlave::_communication->isConnected());


  //necessary as this test does not call finalize
  utils::MasterSlave::_communication = nullptr;
  utils::Parallel::clearGroups();
}

BOOST_AUTO_TEST_CASE(TestFinalize, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "config1.xml";
  if(utils::Parallel::getProcessRank()<=1){
    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 2 );
    slimConfigure(interface, configFilename);
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
    slimConfigure(interface, configFilename);
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

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator()); //needed since this test uses PETSc
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    slimConfigure(interface, configFilename);
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
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverTwo", 0, 1 );
    slimConfigure(interface, configFilename);
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

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3 );
    slimConfigure(interface, configFilename);
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
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverTwo", 0, 1 );
    slimConfigure(interface, configFilename);
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

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "Ateles" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    SolverInterface interface ( "Ateles", utils::Parallel::getProcessRank(), 3 );
    slimConfigure(interface, configFilename);
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
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    SolverInterface interface ( "FASTEST", 0, 1 );
    slimConfigure(interface, configFilename);
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
  BOOST_TEST(utils::Parallel::getCommunicatorSize() == size);
  utils::Parallel::clearGroups();

  for(int k=0; k<numberOfTests; k++){
    reset();

    SolverInterface interface ( solverName, rank, size );
    slimConfigure(interface, configs[k]);
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
    slimConfigure(precice, _pathToTests + fileName);
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

/// This testcase is based on a bug documented in issue #371
BOOST_AUTO_TEST_CASE(NearestProjectionRePartitioning, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "np-repartitioning.xml";

  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "FluidSolver" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 3);
    utils::Parallel::clearGroups();
    SolverInterface interface ( "FluidSolver", utils::Parallel::getProcessRank(), 3 );
    slimConfigure(interface, configFilename);

    if(utils::Parallel::getProcessRank()==1) {

      const int meshID = interface.getMeshID("CellCenters");
      const int dimensions = 3;
      BOOST_TEST(interface.getDimensions()==dimensions);

      const int numberOfVertices = 65;
      const double yCoord = 0.0;
      const double zCoord = 0.005;
      const std::vector<double> positions{
                                         0.00124795, yCoord, zCoord,
                                         0.00375646, yCoord, zCoord,
                                         0.00629033, yCoord, zCoord,
                                         0.00884982, yCoord, zCoord,
                                         0.0114352, yCoord, zCoord,
                                         0.0140467, yCoord, zCoord,
                                         0.0166846, yCoord, zCoord,
                                         0.0193492, yCoord, zCoord,
                                         0.0220407, yCoord, zCoord,
                                         0.0247594, yCoord, zCoord,
                                         0.0275056, yCoord, zCoord,
                                         0.0302796, yCoord, zCoord,
                                         0.0330816, yCoord, zCoord,
                                         0.0359119, yCoord, zCoord,
                                         0.0387709, yCoord, zCoord,
                                         0.0416588, yCoord, zCoord,
                                         0.0445758, yCoord, zCoord,
                                         0.0475224, yCoord, zCoord,
                                         0.0504987, yCoord, zCoord,
                                         0.0535051, yCoord, zCoord,
                                         0.0565419, yCoord, zCoord,
                                         0.0596095, yCoord, zCoord,
                                         0.062708, yCoord, zCoord,
                                         0.0658378, yCoord, zCoord,
                                         0.0689993, yCoord, zCoord,
                                         0.0721928, yCoord, zCoord,
                                         0.0754186, yCoord, zCoord,
                                         0.0786769, yCoord, zCoord,
                                         0.0819682, yCoord, zCoord,
                                         0.0852928, yCoord, zCoord,
                                         0.088651, yCoord, zCoord,
                                         0.0920431, yCoord, zCoord,
                                         0.0954695, yCoord, zCoord,
                                         0.0989306, yCoord, zCoord,
                                         0.102427, yCoord, zCoord,
                                         0.105958, yCoord, zCoord,
                                         0.109525, yCoord, zCoord,
                                         0.113128, yCoord, zCoord,
                                         0.116768, yCoord, zCoord,
                                         0.120444, yCoord, zCoord,
                                         0.124158, yCoord, zCoord,
                                         0.127909, yCoord, zCoord,
                                         0.131698, yCoord, zCoord,
                                         0.135525, yCoord, zCoord,
                                         0.139391, yCoord, zCoord,
                                         0.143296, yCoord, zCoord,
                                         0.147241, yCoord, zCoord,
                                         0.151226, yCoord, zCoord,
                                         0.15525, yCoord, zCoord,
                                         0.159316, yCoord, zCoord,
                                         0.163422, yCoord, zCoord,
                                         0.16757, yCoord, zCoord,
                                         0.17176, yCoord, zCoord,
                                         0.175993, yCoord, zCoord,
                                         0.180268, yCoord, zCoord,
                                         0.184586, yCoord, zCoord,
                                         0.188948, yCoord, zCoord,
                                         0.193354, yCoord, zCoord,
                                         0.197805, yCoord, zCoord,
                                         0.202301, yCoord, zCoord,
                                         0.206842, yCoord, zCoord,
                                         0.211429, yCoord, zCoord,
                                         0.216062, yCoord, zCoord,
                                         0.220742, yCoord, zCoord,
                                         0.22547, yCoord, zCoord};
      BOOST_TEST(numberOfVertices*dimensions == positions.size());
      std::vector<int> vertexIDs(numberOfVertices);
      interface.setMeshVertices(meshID, numberOfVertices, positions.data(), vertexIDs.data());
      interface.initialize();
      BOOST_TEST(impl(interface).getMeshHandle("Nodes").triangles().size()==15);
      interface.advance(1.0);
      interface.finalize();
    } else {
        interface.initialize();
        interface.advance(1.0);
        interface.finalize();
    }
  }
  else {
    utils::Parallel::splitCommunicator( "SolidSolver" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();
    SolverInterface interface ( "SolidSolver", 0, 1 );
    slimConfigure(interface, configFilename);
    const int meshID = interface.getMeshID("Nodes");
    const int dimensions = 3;
    BOOST_TEST(interface.getDimensions()==dimensions);
    const int numberOfVertices = 34;
    const double yCoord = 0.0;
    const double zCoord1 = 0.0;
    const double zCoord2 = 0.01;
    const std::vector<double> positions{
                                      0.0, yCoord, zCoord2,
                                      0.0, yCoord, zCoord1,
                                      0.03125, yCoord, zCoord2,
                                      0.03125, yCoord, zCoord1,
                                      0.0625, yCoord, zCoord2,
                                      0.0625, yCoord, zCoord1,
                                      0.09375, yCoord, zCoord2,
                                      0.09375, yCoord, zCoord1,
                                      0.125, yCoord, zCoord2,
                                      0.125, yCoord, zCoord1,
                                      0.15625, yCoord, zCoord2,
                                      0.15625, yCoord, zCoord1,
                                      0.1875, yCoord, zCoord2,
                                      0.1875, yCoord, zCoord1,
                                      0.21875, yCoord, zCoord2,
                                      0.21875, yCoord, zCoord1,
                                      0.25, yCoord, zCoord2,
                                      0.25, yCoord, zCoord1,
                                      0.28125, yCoord, zCoord2,
                                      0.28125, yCoord, zCoord1,
                                      0.3125, yCoord, zCoord2,
                                      0.3125, yCoord, zCoord1,
                                      0.34375, yCoord, zCoord2,
                                      0.34375, yCoord, zCoord1,
                                      0.375, yCoord, zCoord2,
                                      0.375, yCoord, zCoord1,
                                      0.40625, yCoord, zCoord2,
                                      0.40625, yCoord, zCoord1,
                                      0.4375, yCoord, zCoord2,
                                      0.4375, yCoord, zCoord1,
                                      0.46875, yCoord, zCoord2,
                                      0.46875, yCoord, zCoord1,
                                      0.5, yCoord, zCoord2,
                                      0.5, yCoord, zCoord1};
    BOOST_TEST(numberOfVertices*dimensions == positions.size());
    std::vector<int> vertexIDs(numberOfVertices);
    interface.setMeshVertices(meshID, numberOfVertices, positions.data(), vertexIDs.data());

    const int numberOfCells = numberOfVertices / 2 - 1;
    const int numberOfEdges = numberOfCells * 4 + 1;
    std::vector<int> edgeIDs(numberOfEdges);

    for( int i = 0; i < numberOfCells; i++){
      edgeIDs.at(4*i)   = interface.setMeshEdge(meshID, vertexIDs.at(i*2),   vertexIDs.at(i*2+1)); //left
      edgeIDs.at(4*i+1) = interface.setMeshEdge(meshID, vertexIDs.at(i*2),   vertexIDs.at(i*2+2)); //top
      edgeIDs.at(4*i+2) = interface.setMeshEdge(meshID, vertexIDs.at(i*2+1), vertexIDs.at(i*2+3)); //bottom
      edgeIDs.at(4*i+3) = interface.setMeshEdge(meshID, vertexIDs.at(i*2),   vertexIDs.at(i*2+3)); //diagonal
    }
    edgeIDs.at(numberOfEdges-1) = interface.setMeshEdge(meshID, vertexIDs.at(numberOfVertices-2), vertexIDs.at(numberOfVertices-1)); //very right

    for( int i = 0; i < numberOfCells; i++){
      interface.setMeshTriangle(meshID, edgeIDs.at(4*i),   edgeIDs.at(4*i+3), edgeIDs.at(4*i+2)); //left-diag-bottom
      interface.setMeshTriangle(meshID, edgeIDs.at(4*i+1), edgeIDs.at(4*i+3), edgeIDs.at(4*i+4)); //top-diag-right
    }

    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(MasterSockets, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "master-sockets.xml";
  std::string myName, myMeshName;
  int myRank, mySize;
  if(utils::Parallel::getProcessRank()<=2){
    myName = "ParallelSolver";
    myRank = utils::Parallel::getProcessRank();
    mySize = 3;
    myMeshName = "ParallelMesh";
  }
  else{
    myName = "SerialSolver";
    myRank = 0;
    mySize = 1;
    myMeshName = "SerialMesh";
  }
  SolverInterface interface(myName, myRank, mySize);
  slimConfigure(interface, configFilename);
  int meshID = interface.getMeshID(myMeshName);
  double position[2] = {0, 0};
  interface.setMeshVertex(meshID, position);
  interface.initialize();
  interface.advance(1.0);
  interface.finalize();
}

// Tests SolverInterface() with a user-defined MPI communicator.
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicator, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "userDefinedMPICommunicator.xml";
  
  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    MPI_Comm myComm = utils::Parallel::getLocalCommunicator();
    int myCommSize;
    MPI_Comm_size(myComm, &myCommSize);
    BOOST_TEST(myCommSize == 3);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3, &myComm );
    slimConfigure(interface, configFilename);
    int meshID = interface.getMeshID("MeshOne");

    int vertexIDs[2];
    double xCoord = utils::Parallel::getProcessRank() * 0.4;
    double positions[4] = {xCoord,0.0,xCoord+0.2,0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "SolverTwo" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverTwo", 0, 1 );
    slimConfigure(interface, configFilename);
    int meshID = interface.getMeshID("MeshTwo");
    int vertexIDs[6];
    double positions[12] = {0.0,0.0,0.2,0.0,0.4,0.0,0.6,0.0,0.8,0.0,1.0,0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}

#ifndef PRECICE_NO_PETSC
// Tests SolverInterface() with a user-defined MPI communicator.
// Since PETSc also uses MPI, we use petrbf mapping here.
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicatorPetRBF, * testing::OnSize(4))
{
  std::string configFilename = _pathToTests + "userDefinedMPICommunicatorPetRBF.xml";
  config::Configuration config;
  
  if(utils::Parallel::getProcessRank()<=2){
    utils::Parallel::splitCommunicator( "SolverOne" );
    MPI_Comm myComm = utils::Parallel::getLocalCommunicator();
    int myCommSize;
    MPI_Comm_size(myComm, &myCommSize);
    BOOST_TEST(myCommSize == 3);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverOne", utils::Parallel::getProcessRank(), 3, &myComm );
    slimConfigure(interface, configFilename);
    int meshID = interface.getMeshID("MeshOne");

    int vertexIDs[2];
    double xCoord = utils::Parallel::getProcessRank() * 0.4;
    double positions[4] = {xCoord,0.0,xCoord+0.2,0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
  else {
    utils::Parallel::splitCommunicator( "SolverTwo" );
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
    BOOST_TEST(utils::Parallel::getCommunicatorSize() == 1);
    utils::Parallel::clearGroups();

    SolverInterface interface ( "SolverTwo", 0, 1 );
    slimConfigure(interface, configFilename);
    int meshID = interface.getMeshID("MeshTwo");
    int vertexIDs[6];
    double positions[12] = {0.0,0.0,0.2,0.0,0.4,0.0,0.6,0.0,0.8,0.0,1.0,0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}
#endif // PRECICE_NO_PETSC

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

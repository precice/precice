#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "partition/ProvidedBoundingBox.hpp"
#include "partition/ReceivedBoundingBox.hpp"
#include "partition/SharedPointer.hpp"
#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "utils/MasterSlave.hpp"
#include "com/CommunicationFactory.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedBoundingBoxTests)

void setupParallelEnvironment(m2n::PtrM2N m2n){
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom =
        com::PtrCommunication(new com::SocketCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //Master-A
    utils::Parallel::splitCommunicator( "Fluid" );
    m2n->acceptMasterConnection ( "Fluid", "SolidMaster");
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master-B
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "Fluid", "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1-B
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2-B
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 2;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }

  if(utils::Parallel::getProcessRank() == 1){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "SolidSlaves", utils::Parallel::getProcessRank());
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 0 , 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 1 , 2 );
  }
}


void tearDownParallelEnvironment(){
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  //utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}


void createNastinMesh2D(mesh::PtrMesh pNastinMesh)
{
  int dimensions = 2;
  assertion(pNastinMesh.use_count() > 0);
  assertion(pNastinMesh->getDimensions() == dimensions);

  if (utils::Parallel::getProcessRank() == 1) {

    Eigen::VectorXd position(dimensions);
    position << 0.10, 0.10;
    pNastinMesh->createVertex(position);
    position << 0.90, 0.90;
    pNastinMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) {    
    // not at interface
    Eigen::VectorXd position(dimensions);
    position << 1.1, 1.1;
    pNastinMesh->createVertex(position);
    position << 1.9, 1.9;
    pNastinMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 3) {

    Eigen::VectorXd position(dimensions);
    position << 2.1, 2.1;
    pNastinMesh->createVertex(position);
    position << 2.9, 2.9;
    pNastinMesh->createVertex(position);
  }
}


BOOST_AUTO_TEST_CASE(TestCommunicateBoundingBox, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 2;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize()==4); 

  // construct send global boundingbox
  mesh::Mesh::BoundingBoxMap sendGlobalBB;
  mesh::Mesh::BoundingBox initialBB;  
  for (int remoteRank = 0; remoteRank < 3; remoteRank++ )
  {
    for (int i=0; i < dimensions; i++) {
      initialBB.push_back(std::make_pair(3-remoteRank-1,3-remoteRank));
    }    
    sendGlobalBB[remoteRank]= initialBB;
    initialBB.clear();
  }
    
  if (utils::Parallel::getProcessRank() == 0)
  {
    int remoteParComSize = 0;
    std::map<int, std::vector<int>> receivedFeedbackMap;
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));     
    m2n->getMasterCommunication()->send(3, 0);     
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(sendGlobalBB, 0 );
    m2n->getMasterCommunication()->receive(remoteParComSize, 0);
    
    for (int i=0; i < remoteParComSize; i++)
    {
      std::vector<int> initialFeedback;
      initialFeedback.push_back(-1);
      receivedFeedbackMap[i]=initialFeedback;
    }    
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveFeedbackMap(receivedFeedbackMap, 0 );

    // test wether we receive correct freedback map!
    BOOST_TEST(receivedFeedbackMap[0][0] == 2);
    BOOST_TEST(receivedFeedbackMap[1][0] == 1);
    BOOST_TEST(receivedFeedbackMap[2][0] == 0);

  }  
  else
  {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));

    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh);   
    pNastinMesh->computeState();

    double safetyFactor = 0.0;
    
    ReceivedBoundingBox part(pSolidzMesh, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);    
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();    
  }  
  tearDownParallelEnvironment();
}



BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

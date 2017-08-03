#include "testing/Testing.hpp"

#include "partition/ProvidedPartition.hpp"

#include "partition/ReceivedPartition.hpp"

#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/GatherScatterComFactory.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ProvidedPartitionTests)


BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate, * testing::OnRanks({0, 1, 2, 3}))
{
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getRestrictedCommunicator({0,1,2,3}));
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));
  com::PtrCommunication masterSlaveCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    utils::Parallel::splitCommunicator( "Fluid" );
    m2n->acceptMasterConnection ( "Fluid", "SolidMaster");
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "Fluid", "SolidMaster" );
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    utils::Parallel::splitCommunicator( "SolidSlaves");
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "SolidSlaves");
  }

  if(utils::Parallel::getProcessRank() == 1){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "SolidSlaves", 0, 1);
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 0, 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 1, 2 );
  }


  int dimensions = 2;
  bool flipNormals = false;
  Eigen::VectorXd offset = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    bool filterFirst = false;
    double safetyFactor = 0.1;

    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size()==6);
    BOOST_TEST(pSolidzMesh->edges().size()==4);

  }
  else{//SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      utils::MasterSlave::_rank = 0;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = false;
      utils::MasterSlave::_masterMode = true;
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      mesh::Vertex& v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5;
      mesh::Vertex& v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1,v2);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      utils::MasterSlave::_rank = 1;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      utils::MasterSlave::_rank = 2;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
      Eigen::VectorXd position(dimensions);
      position << 0.0, 3.5;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      position << 0.0, 5.5;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 0.0, 7.0;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
      pSolidzMesh->createEdge(v4,v5);
      pSolidzMesh->createEdge(v5,v6);
    }

    bool hasToSend = true;
    ProvidedPartition part(hasToSend);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();
    part.compute();


    if(utils::Parallel::getProcessRank() == 2){//Slave1
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2]==6);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2]==6);
    }


  }


  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

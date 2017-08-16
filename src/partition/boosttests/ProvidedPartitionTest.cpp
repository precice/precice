#ifndef PRECICE_NO_MPI
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

void setupParallelEnvironment(m2n::PtrM2N m2n){
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom =
        com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    utils::Parallel::splitCommunicator( "Fluid" );
    m2n->acceptMasterConnection ( "Fluid", "SolidMaster");
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "Fluid", "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 2;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
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
}

void tearDownParallelEnvironment(){
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}


BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate2D, * testing::OnRanks({0, 1, 2, 3}))
{
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getRestrictedCommunicator({0,1,2,3}));
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 2;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    bool filterFirst = false;
    double safetyFactor = 0.1;

    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size()==6);
    BOOST_TEST(pSolidzMesh->edges().size()==4);

    for(int i=0; i<6; i++){
      BOOST_TEST(pSolidzMesh->vertices()[i].getGlobalIndex()==i);
    }
  }
  else{//SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      mesh::Vertex& v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5;
      mesh::Vertex& v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1,v2);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
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


    if(utils::Parallel::getProcessRank() == 1){//master
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2]==6);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
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

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate3D, * testing::OnRanks({0, 1, 2, 3}))
{
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getRestrictedCommunicator({0,1,2,3}));
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 3;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    bool filterFirst = false;
    double safetyFactor = 0.1;

    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size()==6);
    BOOST_TEST(pSolidzMesh->edges().size()==6);
    BOOST_TEST(pSolidzMesh->triangles().size()==2);

    for(int i=0; i<6; i++){
      BOOST_TEST(pSolidzMesh->vertices()[i].getGlobalIndex()==i);
    }
  }
  else{//SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0, 0.0;
      mesh::Vertex& v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5, 1.0;
      mesh::Vertex& v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1,v2);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 0.0, 3.5, 0.1;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, 0.2;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      position << 0.0, 5.5, 0.8;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 0.0, 7.0, 0.4;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position);
      mesh::Edge& e1 = pSolidzMesh->createEdge(v3,v4);
      mesh::Edge& e2 = pSolidzMesh->createEdge(v4,v5);
      mesh::Edge& e3 = pSolidzMesh->createEdge(v5,v3);
      mesh::Edge& e4 = pSolidzMesh->createEdge(v3,v6);
      mesh::Edge& e5 = pSolidzMesh->createEdge(v6,v5);

      pSolidzMesh->createTriangle(e1, e2, e3);
      pSolidzMesh->createTriangle(e4, e5, e3);
    }

    bool hasToSend = true;
    ProvidedPartition part(hasToSend);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();
    part.compute();


    if(utils::Parallel::getProcessRank() == 1){//master
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1]==2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2]==6);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
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

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

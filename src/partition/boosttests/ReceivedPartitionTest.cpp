#include "testing/Testing.hpp"

#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"

#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedPartitionTests)

//TODO split up in two tests
BOOST_AUTO_TEST_CASE(TestRePartitionProjectionMapping, * testing::OnRanks({0, 1, 2, 3}))
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

  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    utils::Parallel::splitCommunicator( "Solid" );
    m2n->acceptMasterConnection ( "Solid", "FluidMaster");
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master
    utils::Parallel::splitCommunicator( "FluidMaster" );
    m2n->requestMasterConnection ( "Solid", "FluidMaster");
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    utils::Parallel::splitCommunicator( "FluidSlaves");
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "FluidSlaves");
  }

  if(utils::Parallel::getProcessRank() == 1){//Master
    masterSlaveCom->acceptConnection ( "FluidMaster", "FluidSlaves", 0, 1);
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    masterSlaveCom->requestConnection( "FluidMaster", "FluidSlaves", 0, 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "FluidMaster", "FluidSlaves", 1, 2 );
  }

  int dimensions = 2;
  bool flipNormals = false;
  Eigen::VectorXd offset = Eigen::VectorXd::Zero(dimensions);



  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
    mesh::PtrMesh pSolidzMesh1(new mesh::Mesh("SolidzMesh1", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh2(new mesh::Mesh("SolidzMesh2", dimensions, flipNormals));

    Eigen::VectorXd position(dimensions);

    position << 0.0, 0.0;
    mesh::Vertex& v1_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v1_2 = pSolidzMesh2->createVertex(position);
    position << 0.0, 1.95;
    mesh::Vertex& v2_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v2_2 = pSolidzMesh2->createVertex(position);
    position << 0.0, 2.1;
    mesh::Vertex& v3_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v3_2 = pSolidzMesh2->createVertex(position);
    position << 0.0, 4.5;
    mesh::Vertex& v4_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v4_2 = pSolidzMesh2->createVertex(position);
    position << 0.0, 5.95;
    mesh::Vertex& v5_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v5_2 = pSolidzMesh2->createVertex(position);
    position << 0.0, 6.1;
    mesh::Vertex& v6_1 = pSolidzMesh1->createVertex(position);
    mesh::Vertex& v6_2 = pSolidzMesh2->createVertex(position);
    pSolidzMesh1->createEdge(v1_1,v2_1);
    pSolidzMesh1->createEdge(v2_1,v3_1);
    pSolidzMesh1->createEdge(v3_1,v4_1);
    pSolidzMesh1->createEdge(v4_1,v5_1);
    pSolidzMesh1->createEdge(v5_1,v6_1);
    pSolidzMesh2->createEdge(v1_2,v2_2);
    pSolidzMesh2->createEdge(v2_2,v3_2);
    pSolidzMesh2->createEdge(v3_2,v4_2);
    pSolidzMesh2->createEdge(v4_2,v5_2);
    pSolidzMesh2->createEdge(v5_2,v6_2);

    bool hasToSend = true;
    ProvidedPartition part1(hasToSend);
    ProvidedPartition part2(hasToSend);

    part1.setMesh(pSolidzMesh1);
    part1.setm2n(m2n);
    part1.communicate();

    part2.setMesh(pSolidzMesh2);
    part2.setm2n(m2n);
    part2.communicate();
  }
  else{
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh1(new mesh::Mesh("SolidzMesh1", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh2(new mesh::Mesh("SolidzMesh2", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping1 = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions) );
    mapping::PtrMapping boundingToMapping1 = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions) );
    mapping::PtrMapping boundingFromMapping2 = mapping::PtrMapping (
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions) );
    mapping::PtrMapping boundingToMapping2 = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions) );
    boundingFromMapping1->setMeshes(pSolidzMesh1,pNastinMesh);
    boundingToMapping1->setMeshes(pNastinMesh,pSolidzMesh1);
    boundingFromMapping2->setMeshes(pSolidzMesh2,pNastinMesh);
    boundingToMapping2->setMeshes(pNastinMesh,pSolidzMesh2);

    if(utils::Parallel::getProcessRank() == 1){//Master
      utils::MasterSlave::_rank = 0;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = false;
      utils::MasterSlave::_masterMode = true;
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      pNastinMesh->createVertex(position);
      position << 0.0, 2.0;
      pNastinMesh->createVertex(position);
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
      position << 0.0, 4.0;
      pNastinMesh->createVertex(position);
      position << 0.0, 6.0;
      pNastinMesh->createVertex(position);
    }

    pNastinMesh->computeState();

    bool filterFirst = false;
    double safetyFactor = 0.1;

    ReceivedPartition part1(filterFirst, dimensions, safetyFactor);
    part1.setMesh(pSolidzMesh1);
    part1.setm2n(m2n);

    //TODO change filterFirst to true once implemented
    ReceivedPartition part2(filterFirst, dimensions, safetyFactor);
    part2.setMesh(pSolidzMesh2);
    part2.setm2n(m2n);

    part1.setFromMapping(boundingFromMapping1);
    part1.setToMapping(boundingToMapping1);
    part2.setFromMapping(boundingFromMapping2);
    part2.setToMapping(boundingToMapping2);


    part1.communicate();
    part1.compute();
    part2.communicate();
    part2.compute();

    // check if the sending and filtering worked right
    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(pSolidzMesh1->vertices().size()==2);
      BOOST_TEST(pSolidzMesh1->edges().size()==1);
      BOOST_TEST(pSolidzMesh2->vertices().size()==3);
      BOOST_TEST(pSolidzMesh2->edges().size()==2);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      BOOST_TEST(pSolidzMesh1->vertices().size()==0);
      BOOST_TEST(pSolidzMesh1->edges().size()==0);
      BOOST_TEST(pSolidzMesh2->vertices().size()==0);
      BOOST_TEST(pSolidzMesh2->edges().size()==0);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      BOOST_TEST(pSolidzMesh1->vertices().size()==2);
      BOOST_TEST(pSolidzMesh1->edges().size()==1);
      BOOST_TEST(pSolidzMesh2->vertices().size()==3);
      BOOST_TEST(pSolidzMesh2->edges().size()==2);
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

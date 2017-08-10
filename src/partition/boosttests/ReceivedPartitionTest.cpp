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

void setupParallelEnvironment(m2n::PtrM2N m2n){
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    utils::Parallel::splitCommunicator( "Solid" );
    m2n->acceptMasterConnection ( "Solid", "FluidMaster");
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master
    utils::Parallel::splitCommunicator( "FluidMaster" );
    m2n->requestMasterConnection ( "Solid", "FluidMaster");
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    utils::Parallel::splitCommunicator( "FluidSlaves");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "FluidSlaves");
    utils::MasterSlave::_rank = 2;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
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
}

void tearDownParallelEnvironment(){
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}

void createSolidzMesh2D(mesh::PtrMesh pSolidzMesh){
  int dimensions = 2;
  assertion(pSolidzMesh.use_count()>0);
  assertion(pSolidzMesh->getDimensions()==dimensions);
  Eigen::VectorXd position(dimensions);

  position << 0.0, 0.0;
  mesh::Vertex& v1 = pSolidzMesh->createVertex(position);
  position << 0.0, 1.95;
  mesh::Vertex& v2 = pSolidzMesh->createVertex(position);
  position << 0.0, 2.1;
  mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
  position << 0.0, 4.5;
  mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
  position << 0.0, 5.95;
  mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
  position << 0.0, 6.1;
  mesh::Vertex& v6 = pSolidzMesh->createVertex(position);
  pSolidzMesh->createEdge(v1,v2);
  pSolidzMesh->createEdge(v2,v3);
  pSolidzMesh->createEdge(v3,v4);
  pSolidzMesh->createEdge(v4,v5);
  pSolidzMesh->createEdge(v5,v6);
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh){
  int dimensions = 2;
  assertion(pNastinMesh.use_count()>0);
  assertion(pNastinMesh->getDimensions()==dimensions);

  if(utils::Parallel::getProcessRank() == 1){//Master

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 2.0;
    pNastinMesh->createVertex(position);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    // slave1 not at interface
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2

    Eigen::VectorXd position(dimensions);
    position << 0.0, 4.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 6.0;
    pNastinMesh->createVertex(position);
  }
}


void createSolidzMesh3D(mesh::PtrMesh pSolidzMesh){
  int dimensions = 3;
  Eigen::VectorXd position(dimensions);
  assertion(pSolidzMesh.use_count()>0);
  assertion(pSolidzMesh->getDimensions()==dimensions);

  position << 0.0, 0.0, -0.1;
  mesh::Vertex& v1 = pSolidzMesh->createVertex(position);
  position << -1.0, 0.0, 0.0;
  mesh::Vertex& v2 = pSolidzMesh->createVertex(position);
  position << 1.0, 0.0, 0.0;
  mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
  position << 0.0, -1.0, 0.0;
  mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
  position << 0.0, 1.0, 0.0;
  mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
  mesh::Edge& e1 = pSolidzMesh->createEdge(v1,v2);
  mesh::Edge& e2 = pSolidzMesh->createEdge(v2,v4);
  mesh::Edge& e3 = pSolidzMesh->createEdge(v4,v1);
  mesh::Edge& e4 = pSolidzMesh->createEdge(v1,v3);
  mesh::Edge& e5 = pSolidzMesh->createEdge(v3,v5);
  mesh::Edge& e6 = pSolidzMesh->createEdge(v5,v1);
  pSolidzMesh->createTriangle(e1,e2,e3);
  pSolidzMesh->createTriangle(e4,e5,e6);
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh){
  int dimensions = 3;
  assertion(pNastinMesh.use_count()>0);
  assertion(pNastinMesh->getDimensions()==dimensions);

  if(utils::Parallel::getProcessRank() == 1){//Master

    Eigen::VectorXd position(dimensions);
    position << -1.0, -1.0, 0.0;
    pNastinMesh->createVertex(position);
    position << -0.75, -0.75, 0.5;
    pNastinMesh->createVertex(position);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    // slave1 not at interface
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, -1.0;
    pNastinMesh->createVertex(position);
    position << 0.5, 0.5, 0.0;
    pNastinMesh->createVertex(position);
  }
}


BOOST_AUTO_TEST_CASE(RePartitionNNBroadcastFilter2D, * testing::OnRanks({0, 1, 2, 3}))
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
  Eigen::VectorXd offset = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    createSolidzMesh2D(pSolidzMesh);
    bool hasToSend = true;
    ProvidedPartition part(hasToSend);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();
  }
  else{
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions) );
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions) );
    boundingFromMapping->setMeshes(pSolidzMesh,pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh,pSolidzMesh);

    createNastinMesh2D(pNastinMesh);
    pNastinMesh->computeState();

    bool filterFirst = false;
    double safetyFactor = 0.1;

    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(pSolidzMesh->vertices().size()==2);
      BOOST_TEST(pSolidzMesh->edges().size()==1);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      BOOST_TEST(pSolidzMesh->vertices().size()==0);
      BOOST_TEST(pSolidzMesh->edges().size()==0);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      BOOST_TEST(pSolidzMesh->vertices().size()==2);
      BOOST_TEST(pSolidzMesh->edges().size()==1);
    }

  }

  tearDownParallelEnvironment();
}


BOOST_AUTO_TEST_CASE(RePartitionNPPreFilterPostFilter2D, * testing::OnRanks({0, 1, 2, 3}))
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

  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    createSolidzMesh2D(pSolidzMesh);
    bool hasToSend = true;
    ProvidedPartition part(hasToSend);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();
  }
  else{
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping (
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions) );
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions) );
    boundingFromMapping->setMeshes(pSolidzMesh,pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh,pSolidzMesh);

    createNastinMesh2D(pNastinMesh);

    pNastinMesh->computeState();
    bool filterFirst = false;
    double safetyFactor = 0.1;
    //TODO change filterFirst to true once implemented
    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(pSolidzMesh->vertices().size()==3);
      BOOST_TEST(pSolidzMesh->edges().size()==2);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      BOOST_TEST(pSolidzMesh->vertices().size()==0);
      BOOST_TEST(pSolidzMesh->edges().size()==0);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      BOOST_TEST(pSolidzMesh->vertices().size()==3);
      BOOST_TEST(pSolidzMesh->edges().size()==2);
    }

  }
  tearDownParallelEnvironment();
}


BOOST_AUTO_TEST_CASE(RePartitionNPBroadcastFilter3D, * testing::OnRanks({0, 1, 2, 3}))
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

  if (utils::Parallel::getProcessRank() == 0){ //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    createSolidzMesh3D(pSolidzMesh);
    bool hasToSend = true;
    ProvidedPartition part(hasToSend);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.communicate();
  }
  else{
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping (
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions) );
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping (
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions) );
    boundingFromMapping->setMeshes(pSolidzMesh,pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh,pSolidzMesh);

    createNastinMesh3D(pNastinMesh);

    pNastinMesh->computeState();
    bool filterFirst = false;
    double safetyFactor = 20.0; //should not filter out anything here
    ReceivedPartition part(filterFirst, dimensions, safetyFactor);
    part.setMesh(pSolidzMesh);
    part.setm2n(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(pSolidzMesh->vertices().size()==2);
      BOOST_TEST(pSolidzMesh->edges().size()==1);
      BOOST_TEST(pSolidzMesh->triangles().size()==0);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      BOOST_TEST(pSolidzMesh->vertices().size()==0);
      BOOST_TEST(pSolidzMesh->edges().size()==0);
      BOOST_TEST(pSolidzMesh->triangles().size()==0);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      BOOST_TEST(pSolidzMesh->vertices().size()==3);
      BOOST_TEST(pSolidzMesh->edges().size()==3);
      BOOST_TEST(pSolidzMesh->triangles().size()==1);
    }

  }
  tearDownParallelEnvironment();
}


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

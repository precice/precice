#include <m2n/PointToPointCommunication.hpp>
#include <com/MPIDirectCommunication.hpp>
#include <mesh/Mesh.hpp>

#include <mpi.h>

#include <iostream>

using namespace precice;

int
main(int argc, char** argv) {
  std::cout << "Running communication dummy" << std::endl;

  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  MPI_Comm_size(MPI_COMM_WORLD, &utils::MasterSlave::_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &utils::MasterSlave::_rank);

  if (utils::MasterSlave::_size != 5) {
    std::cout << "Please run with 5 mpi processes" << std::endl;
    return 1;
  }

  if (utils::MasterSlave::_rank == 0) {
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode = false;
  } else {
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode = true;
  }

  if (utils::MasterSlave::_masterMode) {
    utils::Parallel::initialize(NULL, NULL, "Master");
  } else {
    assertion(utils::MasterSlave::_slaveMode);
    utils::Parallel::initialize(NULL, NULL, "Slave");
  }

  utils::MasterSlave::_communication =
      com::PtrCommunication(new com::MPIDirectCommunication());

  int rankOffset = 1;

  if (utils::MasterSlave::_masterMode) {
    utils::MasterSlave::_communication->acceptConnection(
        "Master", "Slave", utils::MasterSlave::_rank, 1);
    utils::MasterSlave::_communication->setRankOffset(rankOffset);
  } else {
    assertion(utils::MasterSlave::_slaveMode);
    utils::MasterSlave::_communication->requestConnection(
        "Master",
        "Slave",
        utils::MasterSlave::_rank - rankOffset,
        utils::MasterSlave::_size - rankOffset);
  }

  mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", 1, true));

  if (utils::MasterSlave::_masterMode) {
    pMesh->setGlobalNumberOfVertices(10);

    pMesh->getVertexDistribution()[0].push_back(1);
    pMesh->getVertexDistribution()[0].push_back(4);

    pMesh->getVertexDistribution()[1].push_back(0);
    pMesh->getVertexDistribution()[1].push_back(2);
    pMesh->getVertexDistribution()[1].push_back(3);

    // pMesh->getVertexDistribution()[3].push_back(3);

    pMesh->getVertexDistribution()[2].push_back(5);
    pMesh->getVertexDistribution()[2].push_back(6);

    pMesh->getVertexDistribution()[4].push_back(7);
    pMesh->getVertexDistribution()[4].push_back(8);
    pMesh->getVertexDistribution()[4].push_back(9);
  }

  m2n::PointToPointCommunication c(pMesh);

  c.acceptConnection("B", "A", 0, 1);

  utils::MasterSlave::_communication.reset();

  MPI_Finalize();
  std::cout << "Stop communication dummy" << std::endl;
}

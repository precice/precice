#include <com/MPIDirectCommunication.hpp>
#include <com/MPIPortsCommunicationFactory.hpp>
#include <com/SocketCommunicationFactory.hpp>
#include <m2n/PointToPointCommunication.hpp>
#include <mesh/Mesh.hpp>
#include <utils/MasterSlave.hpp>

#include <mpi.h>

#include <iostream>
#include <vector>

using namespace precice;

using std::cout;
using std::vector;

vector<double>
getData()
{
  int rank = utils::MasterSlave::getRank();

  static double data_0[] = {10.0, 20.0, 40.0, 80.0};
  static double data_1[] = {30.0, 50.0, 60.0, 90.0};
  static double data_2[] = {70.0, 100.0};

  static double *data[] = {data_0, data_1, data_2};
  static int     size[] = {sizeof(data_0) / sizeof(*data_0),
                       sizeof(data_1) / sizeof(*data_1),
                       sizeof(data_2) / sizeof(*data_2)};

  return std::move(vector<double>(data[rank], data[rank] + size[rank]));
}

vector<double>
getExpectedData()
{
  int rank = utils::MasterSlave::getRank();

  static double data_0[] = {10.0 + 2, 20.0 + 1, 40.0 + 2, 80.0 + 5};
  static double data_1[] = {30.0 + 2, 50.0 + 1, 60.0 + 3, 90.0 + 5};
  static double data_2[] = {70.0 + 3, 100.0 + 5};

  static double *data[] = {data_0, data_1, data_2};
  static int     size[] = {sizeof(data_0) / sizeof(*data_0),
                       sizeof(data_1) / sizeof(*data_1),
                       sizeof(data_2) / sizeof(*data_2)};

  return std::move(vector<double>(data[rank], data[rank] + size[rank]));
}

bool validate(vector<double> const &data)
{
  bool valid = true;

  vector<double> expectedData = getExpectedData();

  if (data.size() != expectedData.size())
    return false;

  for (int i = 0; i < data.size(); ++i) {
    valid &= (data[i] == expectedData[i]);
  }

  return valid;
}

int main(int argc, char **argv)
{
  std::cout << "Running communication dummy\n";

  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  MPI_Comm_size(MPI_COMM_WORLD, &utils::MasterSlave::getSize());
  MPI_Comm_rank(MPI_COMM_WORLD, &utils::MasterSlave::getRank());

  if (utils::MasterSlave::getSize() != 3) {
    std::cout << "Please run with 3 mpi processes\n";
    return 1;
  }

  if (utils::MasterSlave::getRank() == 0) {
    utils::MasterSlave::isMaster() = true;
    utils::MasterSlave::isSlave()  = false;
  } else {
    utils::MasterSlave::isMaster() = false;
    utils::MasterSlave::isSlave()  = true;
  }

  if (utils::MasterSlave::isMaster()) {
    utils::Parallel::initializeMPI(NULL, NULL);
    utils::Parallel::splitCommunicator("Master");
  } else {
    assertion(utils::MasterSlave::isSlave());
    utils::Parallel::initializeMPI(NULL, NULL);
    utils::Parallel::splitCommunicator("Slave");
  }

  utils::MasterSlave::_communication =
      com::PtrCommunication(new com::MPIDirectCommunication);

  int rankOffset = 1;

  if (utils::MasterSlave::isMaster()) {
    utils::MasterSlave::_communication->acceptConnection(
        "Master", "Slave", utils::MasterSlave::getRank(), 1);
    utils::MasterSlave::_communication->setRankOffset(rankOffset);
  } else {
    assertion(utils::MasterSlave::isSlave());
    utils::MasterSlave::_communication->requestConnection(
        "Master",
        "Slave",
        utils::MasterSlave::getRank() - rankOffset,
        utils::MasterSlave::getSize() - rankOffset);
  }

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  if (utils::MasterSlave::isMaster()) {
    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(0);
    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(7);

    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5);
    mesh->getVertexDistribution()[1].push_back(8);

    mesh->getVertexDistribution()[2].push_back(6);
    mesh->getVertexDistribution()[2].push_back(9);
  }

  std::vector<com::PtrCommunicationFactory> cfs(
      {com::PtrCommunicationFactory(new com::SocketCommunicationFactory)});

  //std::vector<com::PtrCommunicationFactory> cfs(
  //     {com::PtrCommunicationFactory(new com::SocketCommunicationFactory),
  //      com::PtrCommunicationFactory(new com::MPIPortsCommunicationFactory)});

  for (auto cf : cfs) {
    m2n::PointToPointCommunication c(cf, mesh);

    c.requestConnection("B", "A");

    cout << utils::MasterSlave::getRank() << ": "
         << "Connected!\n";

    std::vector<double> data = getData();

    c.send(data.data(), data.size());

    c.receive(data.data(), data.size());

    if (validate(data))
      cout << utils::MasterSlave::getRank() << ": "
           << "Success!\n";
    else
      cout << utils::MasterSlave::getRank() << ": "
           << "Failure!\n";

    cout << "----------\n";
  }

  utils::MasterSlave::_communication.reset();

  MPI_Finalize();

  std::cout << "Stop communication dummy\n";
}

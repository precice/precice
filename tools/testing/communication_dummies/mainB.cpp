#include <com/MPIDirectCommunication.hpp>
#include <com/MPIPortsCommunicationFactory.hpp>
#include <com/SocketCommunicationFactory.hpp>
#include <m2n/PointToPointCommunication.hpp>
#include <mesh/Mesh.hpp>
#include <utils/IntraComm.hpp>

#include <mpi.h>

#include <iostream>
#include <vector>

using namespace precice;

using std::cout;
using std::rand;
using std::vector;

vector<double>
getData()
{
  int rank = utils::IntraComm::getRank();

  static double  data_0[] = {rand(), rand()};
  static double  data_1[] = {rand(), rand(), rand()};
  static double  data_2[] = {rand(), rand()};
  static double *data_3;
  static double  data_4[] = {rand(), rand(), rand()};

  static double *data[] = {data_0, data_1, data_2, data_3, data_4};
  static int     size[] = {sizeof(data_0) / sizeof(*data_0),
                           sizeof(data_1) / sizeof(*data_1),
                           sizeof(data_2) / sizeof(*data_2),
                           0,
                           sizeof(data_4) / sizeof(*data_4)};

  return std::move(vector<double>(data[rank], data[rank] + size[rank]));
}

vector<double>
getExpectedData()
{
  int rank = utils::IntraComm::getRank();

  static double  data_0[] = {20.0, 50.0};
  static double  data_1[] = {10.0, 30.0, 40.0};
  static double  data_2[] = {60.0, 70.0};
  static double *data_3;
  static double  data_4[] = {80.0, 90.0, 100.0};

  static double *data[] = {data_0, data_1, data_2, data_3, data_4};
  static int     size[] = {sizeof(data_0) / sizeof(*data_0),
                           sizeof(data_1) / sizeof(*data_1),
                           sizeof(data_2) / sizeof(*data_2),
                           0,
                           sizeof(data_4) / sizeof(*data_4)};

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

void process(vector<double> &data)
{
  for (int i = 0; i < data.size(); ++i) {
    data[i] += utils::IntraComm::getRank() + 1;
  }
}

int main(int argc, char **argv)
{
  std::cout << "Running communication dummy\n";

  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  MPI_Comm_size(MPI_COMM_WORLD, &utils::IntraComm::getSize());
  MPI_Comm_rank(MPI_COMM_WORLD, &utils::IntraComm::getRank());

  if (utils::IntraComm::getSize() != 5) {
    std::cout << "Please run with 5 mpi processes\n";
    return 1;
  }

  if (utils::IntraComm::getRank() == 0) {
    utils::IntraComm::isPrimary()   = true;
    utils::IntraComm::isSecondary() = false;
  } else {
    utils::IntraComm::isPrimary()   = false;
    utils::IntraComm::isSecondary() = true;
  }

  if (utils::IntraComm::isPrimary()) {
    utils::Parallel::initializeMPI(NULL, NULL);
    utils::Parallel::splitCommunicator("Primary");
  } else {
    assertion(utils::IntraComm::isSecondary());
    utils::Parallel::initializeMPI(NULL, NULL);
    utils::Parallel::splitCommunicator("Secondary");
  }

  utils::IntraComm::getCommunication() =
      com::PtrCommunication(new com::MPIDirectCommunication);

  int rankOffset = 1;

  if (utils::IntraComm::isPrimary()) {
    utils::IntraComm::getCommunication()->acceptConnection(
        "Primary", "Secondary", utils::IntraComm::getRank(), 1);
    utils::IntraComm::getCommunication()->setRankOffset(rankOffset);
  } else {
    assertion(utils::IntraComm::isSecondary());
    utils::IntraComm::getCommunication()->requestConnection(
        "Primary",
        "Secondary",
        utils::IntraComm::getRank() - rankOffset,
        utils::IntraComm::getSize() - rankOffset);
  }

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  if (utils::IntraComm::isPrimary()) {
    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(4);

    mesh->getVertexDistribution()[1].push_back(0);
    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(3);

    // mesh->getVertexDistribution()[3].push_back(3);

    mesh->getVertexDistribution()[2].push_back(5);
    mesh->getVertexDistribution()[2].push_back(6);

    mesh->getVertexDistribution()[4].push_back(7);
    mesh->getVertexDistribution()[4].push_back(8);
    mesh->getVertexDistribution()[4].push_back(9);
  }

  std::vector<com::PtrCommunicationFactory> cfs(
      {com::PtrCommunicationFactory(new com::SocketCommunicationFactory)});

  // std::vector<com::PtrCommunicationFactory> cfs(
  //      {com::PtrCommunicationFactory(new com::SocketCommunicationFactory),
  //       com::PtrCommunicationFactory(new com::MPIPortsCommunicationFactory)});

  for (auto cf : cfs) {
    m2n::PointToPointCommunication c(cf, mesh);

    c.acceptConnection("B", "A");

    cout << utils::IntraComm::getRank() << ": "
         << "Connected!" << '\n';

    std::vector<double> data = getData();

    c.receive(data.data(), data.size());

    if (validate(data))
      cout << utils::IntraComm::getRank() << ": "
           << "Success!" << '\n';
    else
      cout << utils::IntraComm::getRank() << ": "
           << "Failure!" << '\n';

    process(data);

    c.send(data.data(), data.size());

    cout << "----------" << '\n';
  }

  utils::IntraComm::getCommunication().reset();

  MPI_Finalize();

  std::cout << "Stop communication dummy\n";
}

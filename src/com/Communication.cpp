#include "Communication.hpp"
#include "Request.hpp"
#include <thread>
#include <chrono>
#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/filesystem.hpp>

#include <iostream>

namespace precice
{
namespace com
{
/**
 * @attention This method modifies the input buffer.
 */
void Communication::reduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);
  
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += itemsToSend[i];
    }
  }
}

void Communication::reduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
}

void Communication::reduceSum(int itemToSend, int &itemToReceive)
{
  TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void Communication::reduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);
  
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += itemsToSend[i];
    }
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToReceive, size, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemsToReceive, size, rankMaster + _rankOffset);
}

void Communication::allreduceSum(double itemToSend, double &itemToReceive)
{
  TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemToSend, 1, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(double itemToSend, double &itemsToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive)
{
  TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemToReceive, 1, rankMaster + _rankOffset);
}

void Communication::broadcast(const int *itemsToSend, int size)
{
  TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(int *itemsToReceive, int size, int rankBroadcaster)
{
  TRACE(size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(int itemToSend)
{
  TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(int &itemToReceive, int rankBroadcaster)
{
  TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(const double *itemsToSend, int size)
{
  TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(double *itemsToReceive,
                              int     size,
                              int     rankBroadcaster)
{
  TRACE(size);
  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(double itemToSend)
{
  TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(double &itemToReceive, int rankBroadcaster)
{
  TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(bool itemToSend)
{
  TRACE();
  int item = itemToSend;
  broadcast(item);
}

void Communication::broadcast(bool &itemToReceive, int rankBroadcaster)
{
  TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void Communication::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(const_cast<int*>(v.data()), v.size()); // make it send vector
}

void Communication::broadcast(std::vector<int> &v, int rankBroadcaster)
{
  v.clear();
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.resize(size);
  broadcast(v.data(), size, rankBroadcaster);
}

void Communication::broadcast(std::vector<double> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(v.data(), v.size()); // make it send vector
}

void Communication::broadcast(std::vector<double> &v, int rankBroadcaster)
{
  v.clear();
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.resize(size);
  broadcast(v.data(), size, rankBroadcaster);
}

void Communication::writeConnectionInfo(std::string const & acceptorName,
                                        std::string const & requesterName,
                                        int rank,
                                        std::string addressDirectory,
                                        std::string addressData)
{
  using namespace boost::filesystem;
  int const firstLevelLen = 2;

  boost::uuids::string_generator ns_gen;
  auto ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace
  
  boost::uuids::name_generator gen{ns};
  std::string s = acceptorName + requesterName + std::to_string(rank);
  std::string hash = boost::uuids::to_string(gen(s));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());

  path p = path(addressDirectory) / path(".precice") / path(hash.substr(0, firstLevelLen));
  create_directories(p);

  p /= hash.substr(firstLevelLen);
  std::cout << "Writing to file:   " << p << " for "
            << acceptorName << ", " << requesterName << ", " << rank << ", " << addressDirectory << std::endl;

  std::ofstream ofs(p.string());
  ofs << addressData;
}


std::string Communication::readConnectionInfo(std::string const & acceptorName,
                                              std::string const & requesterName,
                                              int rank,
                                              std::string addressDirectory)
{
  using namespace boost::filesystem;
  int const firstLevelLen = 2;

  boost::uuids::string_generator ns_gen;
  auto ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace
  
  boost::uuids::name_generator gen{ns};
  std::string s = acceptorName + requesterName + std::to_string(rank);
  std::string hash = boost::uuids::to_string(gen(s));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());
    
  // cout << "Hash = " << hash << endl;
  path p = path(addressDirectory) / path(".precice") / path(hash.substr(0, firstLevelLen));
  
  p /= hash.substr(firstLevelLen);
  std::cout << "Reading from file: " << p << " for "
            << acceptorName << ", " << requesterName << ", " << rank << ", " << addressDirectory << std::endl;

  std::string addressData;
  std::ifstream ifs;
  do {
    ifs.open(p.string(), std::ifstream::in);
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  } while (not ifs);

  ifs >> addressData;
  return addressData;
}

void Communication::removeConnectionInfo(std::string const & acceptorName,
                                      std::string const & requesterName,
                                      int rank,
                                      std::string addressDirectory)
{
  using namespace boost::filesystem;
  int const firstLevelLen = 2;

  boost::uuids::string_generator ns_gen;
  auto ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace
  
  boost::uuids::name_generator gen{ns};
  std::string s = acceptorName + requesterName + std::to_string(rank);
  std::string hash = boost::uuids::to_string(gen(s));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());
  
  path p = path(addressDirectory) / path(".precice") / path(hash.substr(0, firstLevelLen));
  
  p /= hash.substr(firstLevelLen);
  std::cout << "Removing file: " << p << " for "
            << acceptorName << ", " << requesterName << ", " << rank << ", " << addressDirectory << std::endl;

  remove(p);
}



} // namespace com
} // namespace precice

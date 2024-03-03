#ifndef PRECICE_NO_MPI

#include <filesystem>
#include <memory>
#include <mpi.h>
#include <ostream>
#include <utility>

#include "com/ConnectionInfoPublisher.hpp"
#include "com/MPISinglePortsCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/MPIResult.hpp"
#include "utils/Parallel.hpp"
#include "utils/String.hpp"
#include "utils/assertion.hpp"

using precice::utils::MPIResult;

namespace precice::com {
MPISinglePortsCommunication::MPISinglePortsCommunication(std::string addressDirectory)
    : _addressDirectory(std::move(addressDirectory))
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

MPISinglePortsCommunication::~MPISinglePortsCommunication()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

size_t MPISinglePortsCommunication::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  if (_global != MPI_COMM_NULL) {
    int size = -1;
    MPI_Comm_remote_size(_global, &size);
    return size;
  } else {
    return _initialCommSize;
  }
}

void MPISinglePortsCommunication::acceptConnection(std::string const &acceptorName,
                                                   std::string const &requesterName,
                                                   std::string const &tag,
                                                   int                acceptorRank,
                                                   int                rankOffset)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  setRankOffset(rankOffset);

  _isAcceptor = true;

  MPIResult res;

  utils::StringMaker<MPI_MAX_PORT_NAME> sm;
  res = MPI_Open_port(MPI_INFO_NULL, sm.data());
  PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());

  _portName = sm.str();

  ConnectionInfoWriter conPub(acceptorName, requesterName, tag, _addressDirectory);
  conPub.write(_portName);

  int peerCurrent = 0;  // current peer to connect to
  int peerCount   = -1; // The total count of peers (initialized in the first iteration)
  do {
    // Connection
    MPI_Comm communicator;
    res = MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_CHECK(res, "MPI_Comm_accept failed with message: {}", res.message());
    PRECICE_DEBUG("Accepted connection at {} for peer {}", _portName, peerCurrent);

    // Exchange information to which rank I am connected and which communicator size on the other side
    int requesterRank = -1;
    MPI_Recv(&requesterRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    int requesterCommunicatorSize = -1;
    MPI_Recv(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    MPI_Send(&acceptorRank, 1, MPI_INT, 0, 42, communicator);

    // Initialize the count of peers to connect to
    if (peerCurrent == 0) {
      peerCount = requesterCommunicatorSize;
    }

    PRECICE_ASSERT(requesterCommunicatorSize > 0,
                   "Requester communicator size is {} which is invalid.", requesterCommunicatorSize);
    PRECICE_ASSERT(requesterCommunicatorSize == peerCount,
                   "Current requester size from rank {} is {} but should be {}", requesterRank, requesterCommunicatorSize, peerCount);
    PRECICE_ASSERT(_direct.count(requesterRank) == 0,
                   "Rank {} has already been connected. Duplicate requests are not allowed.", requesterRank);

    _direct.emplace(requesterRank, communicator);

    PRECICE_ASSERT(peerCount > 0);
  } while (++peerCurrent < peerCount);

  _initialCommSize = peerCount;
  _isConnected     = true;
}

/// requesterCommunicatorSize is not used, since connection is always made on the entire communicator
void MPISinglePortsCommunication::acceptConnectionAsServer(std::string const &acceptorName,
                                                           std::string const &requesterName,
                                                           std::string const &tag,
                                                           int                acceptorRank,
                                                           int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRank, requesterCommunicatorSize);
  PRECICE_ASSERT(not isConnected());

  _isAcceptor = true;

  const Rank rank = utils::Parallel::current()->rank();
  MPIResult  res;

  if (rank == 0) { // only primary rank opens a port
    ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, _addressDirectory);

    utils::StringMaker<MPI_MAX_PORT_NAME> sm;
    res = MPI_Open_port(MPI_INFO_NULL, sm.data());
    PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());
    _portName = sm.str();

    conInfo.write(_portName);
    PRECICE_DEBUG("Accept connection at {}", _portName);

    res = MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, utils::Parallel::current()->comm, &_global);
    PRECICE_CHECK(res, "MPI_Comm_accept failed with message: {}", res.message());
    PRECICE_DEBUG("Accepted connection at {}", _portName);

  } else { // Secondary ranks call simply call accept

    // The port is only used on the root rank
    res = MPI_Comm_accept(nullptr, MPI_INFO_NULL, 0, utils::Parallel::current()->comm, &_global);
    PRECICE_CHECK(res, "MPI_Comm_accept failed with message: {}", res.message());
    PRECICE_DEBUG("Accepted connection");
  }

  _isConnected = true;
}

void MPISinglePortsCommunication::requestConnection(std::string const &acceptorName,
                                                    std::string const &requesterName,
                                                    std::string const &tag,
                                                    int                requesterRank,
                                                    int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());
  _isAcceptor = false;

  ConnectionInfoReader conInfo(acceptorName, requesterName, tag, _addressDirectory);
  _portName = conInfo.read();
  PRECICE_DEBUG("Request connection to {}", _portName);

  MPI_Comm  communicator;
  MPIResult res = MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
  PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());
  PRECICE_DEBUG("Requested connection to {}", _portName);

  // Send the rank of this requester
  MPI_Send(&requesterRank, 1, MPI_INT, 0, 42, communicator);
  // Send the rank of this requesters communicator size
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);
  // Receive the acceptorRank, which should always be 0
  int acceptorRank = -1;
  MPI_Recv(&acceptorRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  PRECICE_ASSERT(acceptorRank == 0);

  _direct.emplace(acceptorRank, communicator);

  _initialCommSize = requesterCommunicatorSize;
  _isConnected     = true;
}

void MPISinglePortsCommunication::requestConnectionAsClient(std::string const &  acceptorName,
                                                            std::string const &  requesterName,
                                                            std::string const &  tag,
                                                            std::set<int> const &acceptorRanks,
                                                            int                  requesterRank)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  _isAcceptor = false;

  ConnectionInfoReader conInfo(acceptorName, requesterName, tag, _addressDirectory);
  _portName = conInfo.read();
  PRECICE_DEBUG("Request connection to {}", _portName);

  MPIResult res = MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                                   utils::Parallel::current()->comm, &_global);
  PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());
  PRECICE_DEBUG("Requested connection to {}", _portName);

  _isConnected = true;
}

void MPISinglePortsCommunication::closeConnection()
{
  PRECICE_TRACE(_direct.size());

  if (not isConnected())
    return;
  MPIResult res;

  for (auto &kv : _direct) {
    res = MPI_Comm_disconnect(&kv.second);
    if (!res) {
      PRECICE_WARN("MPI_Open_port failed with message: {}", res.message());
    }
  }
  _direct.clear();
  if (_global != MPI_COMM_NULL) {
    res = MPI_Comm_disconnect(&_global);
    if (!res) {
      PRECICE_WARN("MPI_Open_port failed with message: {}", res.message());
    }
  }

  PRECICE_DEBUG("Disconnected");

  if (_isAcceptor and utils::IntraComm::getRank() == 0) {
    res = MPI_Close_port(const_cast<char *>(_portName.c_str()));
    if (!res) {
      PRECICE_WARN("MPI_Open_port failed with message: {}", res.message());
    }
    _portName.clear();
    PRECICE_DEBUG("Port closed");
  }

  _initialCommSize = -1;
  _isConnected     = false;
}

MPI_Comm &MPISinglePortsCommunication::communicator(Rank rank)
{
  if (_global != MPI_COMM_NULL) {
    // Always prefer the global communicator
    return _global;
  } else {
    // Use a direct communication if required
    return _direct.at(rank);
  }
}

int MPISinglePortsCommunication::rank(Rank rank)
{
  if (_global != MPI_COMM_NULL) {
    // Always prefer the global communicator
    return rank;
  } else {
    // Use a direct communication if required.
    // In this case the other rank is always 0.
    return 0;
  }
}

void MPISinglePortsCommunication::prepareEstablishment(std::string const &acceptorName,
                                                       std::string const &requesterName)
{
  using namespace std::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Creating connection exchange directory {}", dir.generic_string());
  try {
    create_directories(dir);
  } catch (const std::filesystem::filesystem_error &e) {
    PRECICE_WARN("Creating directory for connection info failed with: {}", e.what());
  }
}

void MPISinglePortsCommunication::cleanupEstablishment(std::string const &acceptorName,
                                                       std::string const &requesterName)
{
  using namespace std::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Removing connection exchange directory {}", dir.generic_string());
  try {
    remove_all(dir);
  } catch (const std::filesystem::filesystem_error &e) {
    PRECICE_WARN("Cleaning up connection info failed with: {}", e.what());
  }
}

} // namespace precice::com

#endif // not PRECICE_NO_MPI

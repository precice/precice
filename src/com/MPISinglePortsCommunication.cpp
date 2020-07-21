#ifndef PRECICE_NO_MPI

#include "MPISinglePortsCommunication.hpp"
#include <boost/filesystem.hpp>
#include <memory>
#include <mpi.h>
#include <ostream>
#include <utility>
#include "ConnectionInfoPublisher.hpp"
#include "logging/LogMacros.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
MPISinglePortsCommunication::MPISinglePortsCommunication(std::string const &addressDirectory)
    : _addressDirectory(addressDirectory)
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

  MPI_Open_port(MPI_INFO_NULL, const_cast<char *>(_portName.data()));

  ConnectionInfoWriter conPub(acceptorName, requesterName, tag, _addressDirectory);
  conPub.write(_portName);

  int peerCurrent = 0;  // current peer to connect to
  int peerCount   = -1; // The total count of peers (initialized in the first iteration)
  do {
    // Connection
    MPI_Comm communicator;
    MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_DEBUG("Accepted connection at " << _portName << " for peer " << peerCurrent);

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
                   "Requester communicator size is " << requesterCommunicatorSize << " which is invalid.");
    PRECICE_ASSERT(requesterCommunicatorSize == peerCount,
                   "Current requester size from rank " << requesterRank << " is " << requesterCommunicatorSize << " but should be " << peerCount);
    PRECICE_ASSERT(_direct.count(requesterRank) == 0,
                   "Rank " << requesterRank << " has already been connected. Duplicate requests are not allowed.");

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

  const int rank = utils::Parallel::current()->rank();

  if (rank == 0) { // only master opens a port
    ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, _addressDirectory);

    _portName.reserve(MPI_MAX_PORT_NAME);
    MPI_Open_port(MPI_INFO_NULL, const_cast<char *>(_portName.data()));

    conInfo.write(_portName);
    PRECICE_DEBUG("Accept connection at " << _portName);

    MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, utils::Parallel::current()->comm, &_global);
    PRECICE_DEBUG("Accepted connection at " << _portName);

  } else { // Slaves call simply call accept

    // The port is only used on the root rank
    MPI_Comm_accept(nullptr, MPI_INFO_NULL, 0, utils::Parallel::current()->comm, &_global);
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
  PRECICE_DEBUG("Request connection to " << _portName);

  MPI_Comm communicator;
  MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
  PRECICE_DEBUG("Requested connection to " << _portName);

  // Send the rank of this requester
  MPI_Send(&requesterRank, 1, MPI_INT, 0, 42, communicator);
  // Send the rank of this requesters communicator size
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);
  // Recevie the acceptorRank, which should always be 0
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
  PRECICE_DEBUG("Request connection to " << _portName);

  MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                   utils::Parallel::current()->comm, &_global);
  PRECICE_DEBUG("Requested connection to " << _portName);

  _isConnected = true;
}

void MPISinglePortsCommunication::closeConnection()
{
  PRECICE_TRACE(_direct.size());

  if (not isConnected())
    return;

  for (auto &kv : _direct) {
    MPI_Comm_disconnect(&kv.second);
  }
  _direct.clear();
  if (_global != MPI_COMM_NULL) {
    MPI_Comm_disconnect(&_global);
  }

  PRECICE_DEBUG("Disconnected");

  if (_isAcceptor and utils::MasterSlave::getRank() == 0) {
    MPI_Close_port(const_cast<char *>(_portName.c_str()));
    _portName.clear();
    PRECICE_DEBUG("Port closed");
  }

  _initialCommSize = -1;
  _isConnected     = false;
}

MPI_Comm &MPISinglePortsCommunication::communicator(int rank)
{
  if (_global != MPI_COMM_NULL) {
    // Always prefer the global communicator
    return _global;
  } else {
    // Use a direct communication if required
    return _direct.at(rank);
  }
}

int MPISinglePortsCommunication::rank(int rank)
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
  using namespace boost::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Creating connection exchange directory " << dir);
  try {
    create_directories(dir);
  } catch (const boost::filesystem::filesystem_error &e) {
    PRECICE_WARN("Creating directory for connection info failed with: " << e.what());
  }
}

void MPISinglePortsCommunication::cleanupEstablishment(std::string const &acceptorName,
                                                       std::string const &requesterName)
{
  using namespace boost::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Removing connection exchange directory " << dir);
  try {
    remove_all(dir);
  } catch (const boost::filesystem::filesystem_error &e) {
    PRECICE_WARN("Cleaning up connection info failed with: " << e.what());
  }
}

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI

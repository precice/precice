#ifndef PRECICE_NO_MPI

#include "MPISinglePortsCommunication.hpp"
#include <boost/filesystem.hpp>
#include <memory>
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
  int size = -1;
  MPI_Comm_remote_size(_communicators[0], &size);
  return size;
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

  size_t peerCurrent               = 0; // current peer to connect to
  size_t peerCount                 = 0; // The total count of peers (initialized in the first iteration)
  size_t requesterCommunicatorSize = 0;

  do {
    // Connection
    MPI_Comm communicator;
    MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_DEBUG("Accepted connection at " << _portName << " for peer " << peerCurrent);

    int requesterRank = -1;
    // Exchange information to which rank I am connected and which communicator size on the other side
    MPI_Recv(&requesterRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    MPI_Recv(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    MPI_Send(&acceptorRank, 1, MPI_INT, 0, 42, communicator);

    // Initialize the count of peers to connect to
    if (peerCurrent == 0) {
      peerCount = requesterCommunicatorSize;
    }

    PRECICE_ASSERT(requesterCommunicatorSize > 0,
                  "Requester communicator size is " << requesterCommunicatorSize << " which is invalid.");
    PRECICE_ASSERT(requesterCommunicatorSize == peerCount,
                  "Current requester size from rank " << requesterRank << " is " << requesterCommunicatorSize<< " but should be " << peerCount);
    PRECICE_ASSERT(_communicators.count(requesterRank) == 0,
                  "Rank " << requesterRank << " has already been connected. Duplicate requests are not allowed.");

    _communicators[requesterRank] = communicator;

  } while (++peerCurrent < requesterCommunicatorSize);

  _isConnected = true;
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

  ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, _addressDirectory);

  _isAcceptor = true;

  if (utils::MasterSlave::getRank() == 0) { // only master opens a port
    MPI_Open_port(MPI_INFO_NULL, const_cast<char *>(_portName.data()));
    conInfo.write(_portName);
    PRECICE_DEBUG("Accept connection at " << _portName);
  }

  MPI_Comm communicator;
  MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                  utils::Parallel::current()->comm, &communicator);
  PRECICE_DEBUG("Accepted connection at " << _portName);
  _communicators[0] = communicator; // all comms are the same

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

  _isConnected = true;

  int acceptorRank;
  MPI_Send(&requesterRank, 1, MPI_INT, 0, 42, communicator);
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);
  MPI_Recv(&acceptorRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  _communicators[0] = communicator; // should be acceptorRank
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

  MPI_Comm communicator;
  MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                   utils::Parallel::current()->comm, &communicator);
  PRECICE_DEBUG("Requested connection to " << _portName);
  _communicators[0] = communicator; // all comms are the same
  _isConnected      = true;
}

void MPISinglePortsCommunication::closeConnection()
{
  PRECICE_TRACE(_communicators.size());

  if (not isConnected())
    return;

  for (auto &communicator : _communicators) {
    MPI_Comm_disconnect(&communicator.second);
  }

  PRECICE_DEBUG("Disconnected");

  if (_isAcceptor and utils::MasterSlave::getRank() == 0) {
    MPI_Close_port(const_cast<char *>(_portName.c_str()));
    PRECICE_DEBUG("Port closed");
  }

  _isConnected = false;
}

MPI_Comm &MPISinglePortsCommunication::communicator(int rank)
{
  return _communicators[0];
}

int MPISinglePortsCommunication::rank(int rank)
{
  return rank;
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

#ifndef PRECICE_NO_MPI

#include "MPISinglePortsCommunication.hpp"
#include "utils/assertion.hpp"
#include "utils/Parallel.hpp"
#include "utils/Publisher.hpp"

using precice::utils::Publisher;
using precice::utils::ScopedPublisher;

namespace precice
{
namespace com
{
MPISinglePortsCommunication::MPISinglePortsCommunication(std::string const &addressDirectory)
    : _addressDirectory(addressDirectory)
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

MPISinglePortsCommunication::~MPISinglePortsCommunication()
{
  TRACE(_isConnected);
  closeConnection();
}

size_t MPISinglePortsCommunication::getRemoteCommunicatorSize()
{
  TRACE();
  assertion(isConnected());
  return _communicators.size();
}

void MPISinglePortsCommunication::acceptConnection(std::string const &acceptorName,
                                             std::string const &requesterName,
                                             int                acceptorRank)
{
  TRACE(acceptorName, requesterName);
  assertion(not isConnected());

  _isAcceptor = true;

  MPI_Open_port(MPI_INFO_NULL, const_cast<char *>(_portName.data()));

  std::string addressFileName("." + requesterName + "-" + acceptorName + ".address");
  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);
  ScopedPublisher                        p(addressFileName);
  p.write(_portName);
  DEBUG("Accept connection at " << _portName);

  // Connect the first peer, s.t. we can exchange some information about the other side
  MPI_Comm communicator;
  MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
  DEBUG("Accepted connection at " << _portName);

  int    requesterRank      = -1;
  size_t requesterCommunicatorSize = 0;

  // Exchange information to which rank I am connected and which communicator size on the other side
  MPI_Recv(&requesterRank,             1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  MPI_Recv(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  MPI_Send(&acceptorRank,              1, MPI_INT, 0, 42, communicator);
  
  CHECK(requesterCommunicatorSize > 0, "Requester communicator size has to be > 0!");
  _communicators[requesterRank] = communicator;
  
  // Connect all other peers
  for (size_t i = 1; i < requesterCommunicatorSize; ++i) {
    MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    DEBUG("Accepted connection at " << _portName);

    MPI_Recv(&requesterRank,             1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    MPI_Recv(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    MPI_Send(&acceptorRank,              1, MPI_INT, 0, 42, communicator);

    CHECK(requesterCommunicatorSize == _communicators.size(),
          "Requester communicator sizes are inconsistent!");
    CHECK(_communicators.count(requesterRank) == 0,
          "Duplicate request to connect by same rank (" << requesterRank << ")!");

    _communicators[requesterRank] = communicator;
  }
  
  _isConnected = true;
}

void MPISinglePortsCommunication::acceptConnectionAsServer(
    std::string const &acceptorName,
    std::string const &requesterName,
    int                acceptorRank,
    int                requesterCommunicatorSize)
{
  TRACE(acceptorName, requesterName, acceptorRank, requesterCommunicatorSize);
  CHECK(requesterCommunicatorSize > 0, "Requester communicator size has to be > 0!");
  assertion(not isConnected());

  _isAcceptor = true;

  std::string addressFileName("." + requesterName + "-" + acceptorName + ".address");
  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);
  ScopedPublisher                        p(addressFileName);
  if (utils::MasterSlave::_rank == 0) { // only master opens a port
    MPI_Open_port(MPI_INFO_NULL, const_cast<char *>(_portName.data()));
    p.write(_portName);
    DEBUG("Accept connection at " << _portName);
  }
  
  MPI_Comm communicator;
  MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                  utils::Parallel::getGlobalCommunicator(), &communicator);
  DEBUG("Accepted connection at " << _portName);
  _communicators[0] = communicator; // all comms are the same
  _isConnected = true;
}

void MPISinglePortsCommunication::requestConnection(std::string const &acceptorName,
                                              std::string const &requesterName,
                                              int                requesterRank,
                                              int                requesterCommunicatorSize)
{
  TRACE(acceptorName, requesterName);
  assertion(not isConnected());
  _isAcceptor = false;

  std::string addressFileName("." + requesterName + "-" + acceptorName + ".address");
  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);
  Publisher p(addressFileName);
  _portName = p.read();
  DEBUG("Request connection to " << _portName);

  MPI_Comm communicator;
  MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
  DEBUG("Requested connection to " << _portName);

  _isConnected = true;

  int acceptorRank;
  MPI_Send(&requesterRank,             1, MPI_INT, 0, 42, communicator);
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);
  MPI_Recv(&acceptorRank,              1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  _communicators[0] = communicator; // should be acceptorRank
}

void MPISinglePortsCommunication::requestConnectionAsClient(std::string      const &acceptorName,
                                                            std::string      const &requesterName,
                                                            std::set<int>    const &acceptorRanks,
                                                            int                     requesterRank)
                                                      
{
  TRACE(acceptorName, requesterName, acceptorRanks, requesterRank);
  assertion(not isConnected());
  
  _isAcceptor = false;

  std::string addressFileName("." + requesterName + "-" + acceptorName + ".address");
  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);
  Publisher p(addressFileName);
  _portName = p.read();
  DEBUG("Request connection to " << _portName);

  MPI_Comm communicator;
  MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0,
                   utils::Parallel::getGlobalCommunicator(), &communicator);
  DEBUG("Requested connection to " << _portName);
  _communicators[0] = communicator; // all comms are the same
  _isConnected = true;
}

void MPISinglePortsCommunication::closeConnection()
{
  TRACE(_communicators.size());

  if (not isConnected())
    return;

  for (auto communicator : _communicators) {
    MPI_Comm_disconnect(&communicator.second);
  }

  DEBUG("Disconnected");

  if (_isAcceptor and utils::MasterSlave::_rank == 0) {
    MPI_Close_port(const_cast<char *>(_portName.c_str()));
    DEBUG("Port closed");
  }

  _isConnected = false;
}

MPI_Comm &MPISinglePortsCommunication::communicator(int rank)
{
  return _communicators[0];
}

int MPISinglePortsCommunication::rank(int rank)
{
  return 0;
}

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI

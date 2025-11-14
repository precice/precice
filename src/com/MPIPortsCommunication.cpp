#ifndef PRECICE_NO_MPI

#include <filesystem>
#include <ostream>
#include <utility>

#include "com/ConnectionInfoPublisher.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/MPIResult.hpp"
#include "utils/String.hpp"
#include "utils/assertion.hpp"

using precice::utils::MPIResult;

namespace precice::com {
MPIPortsCommunication::MPIPortsCommunication(std::string addressDirectory)
    : _addressDirectory(std::move(addressDirectory))
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

MPIPortsCommunication::~MPIPortsCommunication()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

size_t MPIPortsCommunication::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  return _communicators.size();
}

void MPIPortsCommunication::acceptConnection(std::string const &acceptorName,
                                             std::string const &requesterName,
                                             std::string const &tag,
                                             int                acceptorRank,
                                             int                rankOffset)
{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRank);
  PRECICE_ASSERT(not isConnected());

  setRankOffset(rankOffset);

  _isAcceptor = true;

  MPIResult res;

  Event e0("mpi.acceptConnection.openPort");

  utils::StringMaker<MPI_MAX_PORT_NAME> sm;
  res = MPI_Open_port(MPI_INFO_NULL, sm.data());
  PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());
  _portName = sm.str();

  e0.stop();

  ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, _addressDirectory);
  conInfo.write(_portName);
  PRECICE_DEBUG("Accept connection at {}", _portName);

  int peerCount   = -1; // The total count of peers (initialized in the first iteration)
  int peerCurrent = 0;  // Current peer to connect to
  do {
    Event e1("mpi.acceptConnection.peer." + peerCurrent);
    Event e2("accept");

    // Connection
    MPI_Comm communicator;
    res = MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_CHECK(res, "MPI_Comm_accept failed with message: {}", res.message());
    PRECICE_DEBUG("Accepted connection at {} for peer {}", _portName, peerCurrent);

    e2.stop();
    Event e3("receiveRequesterRankAndCommunicatorSize");

    // Which rank is requesting a connection?
    int requesterRank = -1;
    MPI_Recv(&requesterRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    // How big is the communicator of the requester
    int requesterCommunicatorSize = -1;
    MPI_Recv(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);

    e3.stop();
    Event e4("sendAcceptorRank");

    // Send the rank of the acceptor (this rank).
    MPI_Send(&acceptorRank, 1, MPI_INT, 0, 42, communicator);

    e4.stop();

    // Initialize the count of peers to connect to
    if (peerCurrent == 0) {
      peerCount = requesterCommunicatorSize;
    }

    PRECICE_ASSERT(requesterCommunicatorSize > 0,
                   "Requester communicator size is {} which is invalid.", requesterCommunicatorSize);
    PRECICE_ASSERT(requesterCommunicatorSize == peerCount,
                   "Current requester size from rank {} is {} but should be {}", requesterRank, requesterCommunicatorSize, peerCount);
    PRECICE_ASSERT(_communicators.count(requesterRank) == 0,
                   "Rank {} has already been connected. Duplicate requests are not allowed.", requesterRank);

    _communicators.emplace(requesterRank, communicator);

    e1.stop();
  } while (++peerCurrent < peerCount);

  Event e4("mpi.acceptConnection.closePort");

  res = MPI_Close_port(const_cast<char *>(_portName.c_str()));
  PRECICE_CHECK(res, "MPI_Close_port failed with message: {}", res.message());
  _portName.clear();
  PRECICE_DEBUG("Closed Port");

  e4.stop();

  _isConnected = true;
}

void MPIPortsCommunication::acceptConnectionAsServer(std::string const &acceptorName,
                                                     std::string const &requesterName,
                                                     std::string const &tag,
                                                     int                acceptorRank,
                                                     int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRank, requesterCommunicatorSize);
  PRECICE_ASSERT(requesterCommunicatorSize >= 0, "Requester communicator size has to be positive.");
  PRECICE_ASSERT(not isConnected());

  _isAcceptor = true;
  MPIResult res;

  Event e0("mpi.acceptConnectionAsServer.openPort");

  utils::StringMaker<MPI_MAX_PORT_NAME> sm;
  res = MPI_Open_port(MPI_INFO_NULL, sm.data());
  PRECICE_CHECK(res, "MPI_Open_port failed with message: {}", res.message());
  _portName = sm.str();

  e0.stop();

  ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, acceptorRank, _addressDirectory);
  conInfo.write(_portName);
  PRECICE_DEBUG("Accept connection at {}", _portName);

  for (int connection = 0; connection < requesterCommunicatorSize; ++connection) {
    Event e1("mpi.acceptConnectionAsServer.connection." + connection);
    Event e2("accept");

    MPI_Comm communicator;
    res = MPI_Comm_accept(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_CHECK(res, "MPI_Comm_accept failed with message: {}", res.message());
    PRECICE_DEBUG("Accepted connection at {}", _portName);

    e2.stop();
    Event e3("receiveRequesterRank");

    // Receive the rank of requester
    int requesterRank = -1;
    MPI_Recv(&requesterRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
    PRECICE_ASSERT(requesterRank >= 0, "Invalid requester rank!");

    PRECICE_ASSERT(_communicators.count(requesterRank) == 0, "This connection has already been established.");

    _communicators.emplace(requesterRank, communicator);

    e3.stop();
    e1.stop();
  }

  Event e4("mpi.acceptConnectionAsServer.closePort");

  res = MPI_Close_port(const_cast<char *>(_portName.c_str()));
  PRECICE_CHECK(res, "MPI_Close_port failed with message: {}", res.message());
  _portName.clear();
  PRECICE_DEBUG("Closed Port");

  e4.stop();

  _isConnected = true;
}

void MPIPortsCommunication::requestConnection(std::string const &acceptorName,
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

  Event e0("connect");

  MPI_Comm  communicator;
  MPIResult res = MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
  PRECICE_CHECK(res, "MPI_Comm_connect failed with message: {}", res.message());
  PRECICE_DEBUG("Requested connection to {}", _portName);

  _isConnected = true;

  e0.stop();
  Event e1("sendRequesterRankAndCommunicatorSize");

  // Send the rank of the requester (this rank)
  MPI_Send(&requesterRank, 1, MPI_INT, 0, 42, communicator);
  // Send the size of the requester communicator size
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);

  e1.stop();
  Event e2("receiveAcceptorRank");

  // Receive the rank of the acceptor that we connected to.
  int acceptorRank = -1;
  MPI_Recv(&acceptorRank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  // @todo The following assertion should always be the case, however the
  // acceleration package currently violates this in order to create a circular
  // intra Communication.
  //
  // PRECICE_ASSERT(acceptorRank == 0, "The acceptor always has to be 0.");

  e2.stop();

  _communicators.emplace(acceptorRank, communicator);
}

void MPIPortsCommunication::requestConnectionAsClient(std::string const   &acceptorName,
                                                      std::string const   &requesterName,
                                                      std::string const   &tag,
                                                      std::set<int> const &acceptorRanks,
                                                      int                  requesterRank)

{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRanks, requesterRank);
  PRECICE_ASSERT(not isConnected());

  _isAcceptor = false;

  for (int acceptorRank : acceptorRanks) {
    Event e("mpi.requestConnectionAsClient." + acceptorRank);

    ConnectionInfoReader conInfo(acceptorName, requesterName, tag, acceptorRank, _addressDirectory);
    _portName = conInfo.read();
    PRECICE_DEBUG("Request connection to {}", _portName);

    Event e0("connect");

    MPI_Comm  communicator;
    MPIResult res = MPI_Comm_connect(const_cast<char *>(_portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);
    PRECICE_CHECK(res, "MPI_Comm_connect failed with message: {}", res.message());
    PRECICE_DEBUG("Requested connection to {}", _portName);

    e0.stop();
    Event e1("send");

    // Rank 0 is always the peer, because we connected on COMM_SELF
    MPI_Send(&requesterRank, 1, MPI_INT, 0, 42, communicator);

    e1.stop();
    Event e2("storeCommunicator");

    PRECICE_ASSERT(_communicators.count(acceptorRank) == 0, "This connection has already been established.");
    _communicators.emplace(acceptorRank, communicator);

    e2.stop();
    e.stop();
  }
  _isConnected = true;
}

void MPIPortsCommunication::closeConnection()
{
  PRECICE_TRACE(_communicators.size());

  if (not isConnected())
    return;

  for (auto &communicator : _communicators) {
    MPIResult res = MPI_Comm_disconnect(&communicator.second);
    PRECICE_WARN_IF(!res,
                    "MPI_Comm_disconnect failed with message: {}", res.message());
  }
  _communicators.clear();

  PRECICE_DEBUG("Disconnected");

  _isConnected = false;
}

void MPIPortsCommunication::prepareEstablishment(std::string const &acceptorName,
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

void MPIPortsCommunication::cleanupEstablishment(std::string const &acceptorName,
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

MPI_Comm &MPIPortsCommunication::communicator(Rank rank)
{
  PRECICE_TRACE(rank, _communicators.size(), _isAcceptor);
  // Use bounds checking here, because a std::map otherwise creates element
  return _communicators.at(rank);
}

int MPIPortsCommunication::rank(Rank rank)
{
  return 0;
}

} // namespace precice::com

#endif // not PRECICE_NO_MPI

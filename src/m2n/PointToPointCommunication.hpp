#pragma once

#include "DistributedCommunication.hpp"
#include <list>
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice
{
namespace m2n
{
/**
 * @brief Point-to-point communication implementation of DistributedCommunication.
 *
 * Direct communication of local data subsets is performed between processes of
 * coupled participants. The two supported implementations of direct
 * communication are SocketCommunication and MPIPortsCommunication which can be
 * supplied via their corresponding instantiation factories
 * SocketCommunicationFactory and MPIPortsCommunicationFactory.
 *
 * For the detailed implementation documentation refer to PointToPointCommunication.cpp.
 */
class PointToPointCommunication : public DistributedCommunication
{
public:
  struct ScopedSetEventNamePrefix {
    explicit ScopedSetEventNamePrefix(std::string const &prefix);

    ~ScopedSetEventNamePrefix();

  private:
    std::string _prefix;
  };

public:
  static void setEventNamePrefix(std::string const &prefix);

  static std::string const &eventNamePrefix();

public:
  PointToPointCommunication(com::PtrCommunicationFactory communicationFactory,
                            mesh::PtrMesh                mesh);

  virtual ~PointToPointCommunication();

  /// Returns true, if a connection to a remote participant has been established.
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param[in] acceptorName  Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName);

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /**
   * @brief Sends a subset of local double values corresponding to local indices
   *        deduced from the current and remote vertex distributions.
   */
  virtual void send(double *itemsToSend, size_t size, int valueDimension = 1);

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  virtual void receive(double *itemsToReceive,
                       size_t  size,
                       int     valueDimension = 1);

private:

  /// Checks all stored requests for completion and removes associated buffers
  /**
   * @param[in] blocking False means that the function returns, even when there are requests left.
   */  
  void checkBufferedRequests(bool blocking);
  
  logging::Logger _log{"m2n::PointToPointCommunication"};

  static std::string _prefix;

  com::PtrCommunicationFactory _communicationFactory;

  /**
   * @brief Defines mapping between:
   *        1. local (to the current process) remote process rank;
   *        2. global remote process rank;
   *        3. local data indices, which define a subset of local (for process
   *           rank in the current participant) data to be communicated between
   *           the current process rank and the remote process rank;
   *        4. communication object (provides point-to-point communication
   *           routines).
   */
  struct Mapping {
    int                   localRemoteRank;
    int                   globalRemoteRank;
    std::vector<int>      indices;
    com::PtrCommunication communication;
    com::PtrRequest       request;
    size_t                offset;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        mappings (one to service each point-to-point connection).
   */
  std::vector<Mapping> _mappings;

  std::vector<double> _buffer;

  size_t _localIndexCount = 0;

  size_t _totalIndexCount = 0;

  bool _isConnected = false;

  std::list<std::pair<std::shared_ptr<com::Request>,
                      std::shared_ptr<std::vector<double>>>> bufferedRequests;

};
} // namespace m2n
} // namespace precice

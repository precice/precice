// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_
#define PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_

#include "DistributedCommunication.hpp"

#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace m2n {
/**
 * @brief Point-to-point communication implementation of
 *        DistributedCommunication.
 *
 * Direct communication of local data subsets is performed between processes of
 * coupled participants. The two supported implementations of direct
 * communication are SocketCommunication and MPIPortsCommunication which can be
 * supplied via their corresponding instantiation factories
 * SocketCommunicationFactory and MPIPortsCommunicationFactory.
 *
 * For the detailed implementation documentation refer to
 * PointToPointCommunication.cpp.
 */
class PointToPointCommunication : public DistributedCommunication {
public:
  struct ScopedSetEventNamePrefix {
    ScopedSetEventNamePrefix(std::string const& prefix);

    ~ScopedSetEventNamePrefix();

  private:
    std::string _prefix;
  };

public:
  static void setEventNamePrefix(std::string const& prefix);

  static std::string const& eventNamePrefix();

public:
  /**
   * @brief Constructor.
   */
  PointToPointCommunication(
      com::CommunicationFactory::SharedPointer communicationFactory,
      mesh::PtrMesh mesh);

  /**
   * @brief Destructor.
   */
  virtual ~PointToPointCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been
   *        established.
   */
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection(std::string const& nameAcceptor,
                                std::string const& nameRequester);

  /**
   * @brief Requests connection from participant, which has to call
   *        acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester);

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
  virtual void send(double* itemsToSend, int size, int valueDimension = 1);

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  virtual void receive(double* itemsToReceive,
                       int size,
                       int valueDimension = 1);

private:
  static tarch::logging::Log _log;

  static std::string _prefix;

private:
  com::CommunicationFactory::SharedPointer _communicationFactory;

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
    int localRemoteRank;
    int globalRemoteRank;
    std::vector<int> indices;
    com::Communication::SharedPointer communication;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        mappings (one to service each point-to-point connection).
   */
  std::vector<Mapping> _mappings;

  std::vector<com::Request::SharedPointer> _requests;

  std::vector<double> _buffer;

  size_t _localIndexCount;

  size_t _totalIndexCount;

  bool _isConnected;
};
}
} // namespace precice, m2n

#endif /* PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_ */

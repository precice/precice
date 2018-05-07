#ifndef PRECICE_M2N_GATHER_SCATTER_COMMUNICATION_HPP_
#define PRECICE_M2N_GATHER_SCATTER_COMMUNICATION_HPP_

#include "DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace m2n
{

/**
 * @brief Implements DistributedCommunication by using a gathering/scattering methodology.
 * Arrays of data are always gathered and scattered at the master. No direct communication
 * between slaves is used.
 * For more details see m2n/DistributedCommunication.hpp
 */
class GatherScatterCommunication : public DistributedCommunication
{
public:
  GatherScatterCommunication(
      com::PtrCommunication com,
      mesh::PtrMesh         mesh);

  virtual ~GatherScatterCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call requestConnection().
   *
   * If several connections are going in to a server, the server has to call this
   * method, while the clients have to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection(
      const std::string &nameAcceptor,
      const std::string &nameRequester);

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * If several connections are going in to a server, the clients have to call this
   * method, while the server has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection(
      const std::string &nameAcceptor,
      const std::string &nameRequester);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /// Sends an array of double values from all slaves (different for each slave).
  virtual void send(
      double *itemsToSend,
      size_t  size,
      int     valueDimension);

  /// All slaves receive an array of doubles (different for each slave).
  virtual void receive(
      double *itemsToReceive,
      size_t  size,
      int     valueDimension);

private:
  logging::Logger _log{"m2n::GatherScatterCommunication"};

  /// master to master basic communication
  com::PtrCommunication _com;

  /// Global communication is set up or not
  bool _isConnected;
};

} // namespace m2n
} // namespace precice

#endif /* PRECICE_M2N_GATHER_SCATTER_COMMUNICATION_HPP_ */

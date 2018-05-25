#pragma once
#ifndef PRECICE_NO_MPI

#include <vector>
#include "MPICommunication.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace com
{
/**
 * @brief Provides connection methods based on MPI ports (part of MPI 2.0).
 *
 * The two participants to be connected can be run in two process groups started
 * up individually, i.e. not within the same process group.
 */
class MPIPortsCommunication : public MPICommunication
{
public:
  explicit MPIPortsCommunication(std::string const &addressDirectory = ".");

  virtual ~MPIPortsCommunication();

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * @pre A connection to the remote participant has been setup.
   */
  virtual size_t getRemoteCommunicatorSize() override;

  /// See precice::com::Communication::acceptConnection().
  virtual void acceptConnection(std::string const &nameAcceptor,
                                std::string const &nameRequester) override;

  virtual void acceptConnectionAsServer(std::string const &nameAcceptor,
                                        std::string const &nameRequester,
                                        int                acceptorRank,
                                        int                requesterCommunicatorSize) override;
  
  /// See precice::com::Communication::requestConnection().
  virtual void requestConnection(std::string const &nameAcceptor,
                                 std::string const &nameRequester,
                                 int                requesterProcessRank,
                                 int                requesterCommunicatorSize) override;

  virtual int requestConnectionAsClient(std::string const &nameAcceptor,
                                        std::string const &nameRequester,
                                        int                acceptorRank) override;

  /// See precice::com::Communication::closeConnection().
  virtual void closeConnection() override;

private:
  virtual MPI_Comm &communicator(int rank) override;

  virtual int rank(int rank) override;

  logging::Logger _log{"com::MPIPortsCommunication"};

  std::string _addressDirectory;

  std::vector<MPI_Comm> _communicators;

  /// Name of the port used for connection.
  std::string _portName = std::string(MPI_MAX_PORT_NAME, '\0');

  bool _isAcceptor = false;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI

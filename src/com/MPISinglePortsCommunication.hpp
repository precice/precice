#pragma once
#ifndef PRECICE_NO_MPI

#include <map>
#include <mpi.h>
#include <set>
#include <stddef.h>
#include <string>
#include "MPICommunication.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace com {
/**
 * @brief Provides connection methods based on MPI ports (part of MPI 2.0).
 *
 * The two participants to be connected can be run in two process groups started
 * up individually, i.e. not within the same process group.
 *
 * Notes on the implementation:
 *
 * acceptConnection / requestConnection still uses one communicator per connection, created on MPI_COMM_SELF.
 *
 * acceptConnectionAsServer / requestConnectionAsClient just uses communicator[0], so it actually creates just one comm.
 *
 * The reason of that is that the first variant is called in various different ways, mostly just from the accepting/receiving rank.
 * In order to just use a single communicator, accept/request must be called on the entire global communicator.
 * This would require quite some rework of the calling code. The gains would also be minuscule, because accept/request is mainly used for 1:1 connection anyways.
 *
 * The latter, in contrast, is only called from m2n::PointToPointCommunication, due to that the p2p class was also heavily modified. But this connection method is the most important, because it does the m2n heavy lifiting.
 *
 * If we agree, that acceptConnection / requestConnection just does 1:1 connection, we can rewrite and simplifiy the code.
**/
class MPISinglePortsCommunication : public MPICommunication {
public:
  explicit MPISinglePortsCommunication(std::string const &addressDirectory = ".");

  virtual ~MPISinglePortsCommunication();

  virtual size_t getRemoteCommunicatorSize() override;

  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                rankOffset = 0) override;

  virtual void acceptConnectionAsServer(std::string const &acceptorName,
                                        std::string const &requesterName,
                                        std::string const &tag,
                                        int                acceptorRank,
                                        int                requesterCommunicatorSize) override;

  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName,
                                 std::string const &tag,
                                 int                requesterRank,
                                 int                requesterCommunicatorSize) override;

  virtual void requestConnectionAsClient(std::string const &  acceptorName,
                                         std::string const &  requesterName,
                                         std::string const &  tag,
                                         std::set<int> const &acceptorRanks,
                                         int                  requesterRank) override;

  virtual void prepareEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

  virtual void cleanupEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

  virtual void closeConnection() override;

private:
  virtual MPI_Comm &communicator(int rank) override;

  virtual int rank(int rank) override;

  logging::Logger _log{"com::MPISinglePortsCommunication"};

  std::string _addressDirectory;

  /** @brief A map of direct communication channels based on MPI_COMM_SELF on both sides
   *
   * These 1-1 connections are used until the has been a @ref _global communicator established.
   *
   * These direct connections connect MPI_COMM_SELF on both sides.
   * The call to establish such a connection is thus not collective.
   *
   * These connections are required for the Master-Master and Master-Slaves connections.
   */
  std::map<int, MPI_Comm> _direct;

  /** @brief The global inter-communicator that connects all ranks
   *
   * Once established, this is the default communicator for all communication.
   *
   * The call to establish this communicator is collective over both communicators.
   * Thus the call to @ref requestConnectionAsClient() and @ref acceptConnectionAsServer() is the
   * only time window where this global communicator can be established.
   */
  MPI_Comm _global = MPI_COMM_NULL;

  /// The communicator size known from acceptConnection and requestConnection
  int _initialCommSize = -1;

  /// Name of the port used for connection.
  std::string _portName = std::string(MPI_MAX_PORT_NAME, '\0');

  bool _isAcceptor = false;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI

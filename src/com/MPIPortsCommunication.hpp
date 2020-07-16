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
 */
class MPIPortsCommunication : public MPICommunication {
public:
  explicit MPIPortsCommunication(std::string const &addressDirectory = ".");

  virtual ~MPIPortsCommunication();

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

  virtual void closeConnection() override;

  virtual void prepareEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

  virtual void cleanupEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

private:
  virtual MPI_Comm &communicator(int rank) override;

  virtual int rank(int rank) override;

  logging::Logger _log{"com::MPIPortsCommunication"};

  std::string _addressDirectory;

  /// Remote rank -> communicator map
  std::map<int, MPI_Comm> _communicators;

  /// Name of the port used for connection.
  std::string _portName = std::string(MPI_MAX_PORT_NAME, '\0');

  bool _isAcceptor = false;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI

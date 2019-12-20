#pragma once
#include <string>
#include "logging/Logger.hpp"

namespace precice {
namespace com {


class ConnectionInfoPublisher
{
public:
  
  ConnectionInfoPublisher(std::string acceptorName,
                          std::string requesterName,
                          int rank,
                          std::string addressDirectory) noexcept
    :
    acceptorName(std::move(acceptorName)),
    requesterName(std::move(requesterName)),
    rank(rank),
    addressDirectory(std::move(addressDirectory))
  {}


  ConnectionInfoPublisher(std::string acceptorName,
                          std::string requesterName,
                          std::string addressDirectory) noexcept
    :
    acceptorName(std::move(acceptorName)),
    requesterName(std::move(requesterName)),
    addressDirectory(std::move(addressDirectory))
  {}


protected:

  std::string const acceptorName;
  std::string const requesterName;
  int const rank = -1;
  std::string const addressDirectory;

  /// Returns the file name for the connection information.
  /**
   * It has the form addressDirectory/precice-run/first two letters from hash of 
   * (acceptorName, requesterName, rank)/rest of hash.
   */
  std::string getFilename() const;

  mutable logging::Logger _log{"com::ConnectionInfoPublisher"};
};


/// Reads the connection info for the given participant/rank information
class ConnectionInfoReader : public ConnectionInfoPublisher
{
public:
  
  ConnectionInfoReader(std::string acceptorName,
                       std::string requesterName,
                       int rank,
                       std::string addressDirectory) noexcept
    : ConnectionInfoPublisher(acceptorName, requesterName, rank, addressDirectory)
  {}


  ConnectionInfoReader(std::string acceptorName,
                       std::string requesterName,
                       std::string addressDirectory) noexcept
    : ConnectionInfoPublisher(acceptorName, requesterName, addressDirectory)
  {}


  /// Reads the info from the connection info file. Will block, if the the file is not present.
  std::string read() const;
};


/// Writes the connection info for the given participant/rank information.
/**
 * The file is removed, when the object is destroyed.
 */
class ConnectionInfoWriter : public ConnectionInfoPublisher
{
public:
  
  ConnectionInfoWriter(std::string acceptorName,
                       std::string requesterName,
                       int rank,
                       std::string addressDirectory) noexcept
    : ConnectionInfoPublisher(acceptorName, requesterName, rank, addressDirectory)
  {}
  
  
  ConnectionInfoWriter(std::string acceptorName,
                       std::string requesterName,
                       std::string addressDirectory) noexcept
    : ConnectionInfoPublisher(acceptorName, requesterName, addressDirectory)
  {}


  /// Removes the connection info file and the directories ./precice-run/[hash], is empty.
  ~ConnectionInfoWriter();

  /// Write the string info, e.g. IP:port to the connection info file
  /**
   * which is determined by acceptorName, requesterName, rank, addressDirectory
   * set at construction.
   */
  void write(std::string const & info) const;
};


}
}

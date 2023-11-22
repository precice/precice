#pragma once
#include <string>
#include <string_view>
#include <utility>

#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"

namespace precice::com {

namespace impl {
/// Returns the file name for the connection information.
/**
   * It has the form first two letters from hash of
   * (acceptorName, requesterName, mesh, rank)/rest of hash.
   */
std::string hashedFilePath(std::string_view acceptorName, std::string_view requesterName, std::string_view meshName, Rank rank);

/** Returns the local directory which is the root for storing connection information.
   * It has the form addressDirectory/precice-run/acceptorName-requesterName
   */
std::string localDirectory(std::string_view acceptorName, std::string_view requesterName, std::string_view addressDirectory);
} // namespace impl

class ConnectionInfoPublisher {
public:
  ConnectionInfoPublisher(std::string_view acceptorName,
                          std::string_view requesterName,
                          std::string_view tag,
                          int              rank,
                          std::string_view addressDirectory) noexcept
      : acceptorName(acceptorName),
        requesterName(requesterName),
        tag(tag),
        rank(rank),
        addressDirectory(addressDirectory)
  {
  }

  ConnectionInfoPublisher(std::string_view acceptorName,
                          std::string_view requesterName,
                          std::string_view tag,
                          std::string_view addressDirectory) noexcept
      : acceptorName(acceptorName),
        requesterName(requesterName),
        tag(tag),
        addressDirectory(addressDirectory)
  {
  }

protected:
  std::string const acceptorName;
  std::string const requesterName;
  std::string const tag;
  int const         rank = -1;
  std::string const addressDirectory;

  /// Returns the local directory which is used to store the hashed part.
  std::string getLocalDirectory() const;

  /// Returns the full path to the hashed filename
  std::string getFilename() const;

  mutable logging::Logger _log{"com::ConnectionInfoPublisher"};
};

/// Reads the connection info for the given participant/rank information
class ConnectionInfoReader : public ConnectionInfoPublisher {
public:
  ConnectionInfoReader(std::string_view acceptorName,
                       std::string_view requesterName,
                       std::string_view tag,
                       int              rank,
                       std::string_view addressDirectory) noexcept
      : ConnectionInfoPublisher(acceptorName, requesterName, tag, rank, addressDirectory)
  {
  }

  ConnectionInfoReader(std::string_view acceptorName,
                       std::string_view requesterName,
                       std::string_view tag,
                       std::string_view addressDirectory) noexcept
      : ConnectionInfoPublisher(acceptorName, requesterName, tag, addressDirectory)
  {
  }

  /// Reads the info from the connection info file. Will block, if the the file is not present.
  std::string read() const;
};

/// Writes the connection info for the given participant/rank information.
/**
 * The file is removed, when the object is destroyed.
 */
class ConnectionInfoWriter : public ConnectionInfoPublisher {
public:
  ConnectionInfoWriter(std::string_view acceptorName,
                       std::string_view requesterName,
                       std::string_view tag,
                       int              rank,
                       std::string_view addressDirectory) noexcept
      : ConnectionInfoPublisher(acceptorName, requesterName, tag, rank, addressDirectory)
  {
  }

  ConnectionInfoWriter(std::string_view acceptorName,
                       std::string_view requesterName,
                       std::string_view tag,
                       std::string_view addressDirectory) noexcept
      : ConnectionInfoPublisher(acceptorName, requesterName, tag, addressDirectory)
  {
  }

  /// Removes the connection info file and the directories ./precice-run/[hash], is empty.
  ~ConnectionInfoWriter();

  /// Write the string info, e.g. IP:port to the connection info file
  /**
   * which is determined by acceptorName, requesterName, rank, addressDirectory
   * set at construction.
   */
  void write(std::string_view info) const;
};

} // namespace precice::com

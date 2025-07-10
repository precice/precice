#pragma once

#include <string>
#include <string_view>

#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice::m2n {

/// An M2N between participants with a configured direction
class BoundM2N {
public:
  /// Prepare to establish the connection
  void prepareEstablishment();

  /// Connect the Primary Ranks of the M2N
  void connectPrimaryRanks(std::string_view configHash);

  /// Connect the Secondary ranks of the M2N
  void connectSecondaryRanks();

  /// pre-connect the Secondary ranks of the M2N
  void preConnectSecondaryRanks();

  /// Cleanup after having established the connection
  void cleanupEstablishment();

  PtrM2N      m2n;
  std::string localName;
  std::string remoteName;
  bool        isRequesting = false;

private:
  mutable logging::Logger _log{"impl::ParticipantImpl"};

  /** Instructs the Primary rank wait for SecondaryRanks.
   *
   * Performs a collective operation which forces every secondary rank to sync with the Primary.
   *
   * @note this does nothing if the participant is running serially.
   */
  void waitForSecondaryRanks();
};

} // namespace precice::m2n

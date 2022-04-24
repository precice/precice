#pragma once

#include <string>
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice {
namespace m2n {

/// An M2N between participants with a configured direction
class BoundM2N {
public:
  /// Prepare to establish the connection
  void prepareEstablishment();

  /// Connect the Primaries of the M2N
  void connectPrimaries();

  /// Connect the Secondaries of the M2N
  void connectSecondaries();

  /// pre-connect the Secondaries of the M2N
  void preConnectSecondaries();

  /// Cleanup after having established the connection
  void cleanupEstablishment();

  PtrM2N      m2n;
  std::string localName;
  std::string remoteName;
  bool        isRequesting = false;

private:
  mutable logging::Logger _log{"impl::SolverInterfaceImpl"};

  /** Instructs the Primary wait for Secondaries.
   *
   * Performs a collective operation which forces every slave to sync with the Primary.
   * 
   * @note this does nothing if the participant is running serially.
   */
  void waitForSecondaries();
};

} // namespace m2n
} // namespace precice

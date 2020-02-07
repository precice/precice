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

  /// Connect the Masters of the M2N
  void connectMasters();

  /// Connect the Slaves of the M2N
  void connectSlaves();

  /// pre-connect the Slaves of the M2N
  void preConnectSlaves();

  /// Cleanup after having established the connection
  void cleanupEstablishment();

  PtrM2N      m2n;
  std::string localName;
  std::string remoteName;
  bool        isRequesting;

private:
  mutable logging::Logger _log{"impl::SolverInterfaceImpl"};

  /** Instructs the Master wait for Slaves.
   *
   * Performs a collective operation which forces every slave to sync with the Master.
   * 
   * @note this does nothing if the participant is running serially.
   */
  void waitForSlaves();
};

} // namespace m2n
} // namespace precice

#include "com/ErrorHandling.hpp"

namespace precice::com {

void checkErrorCode(boost::system::error_code ec, logging::Logger &_log)
{
  if (ec.value() == boost::system::errc::success) {
    return;
  }
  if (ec.value() == boost::system::errc::connection_reset) {
    PRECICE_ERROR("Connection was reset by another participant, which most likely exited unexpectedly (look there).");
    return;
  }
  PRECICE_ERROR("Receiving data from another participant failed with a system error: {}.", ec.message());
}
} // namespace precice::com

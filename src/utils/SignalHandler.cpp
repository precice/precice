#include "SignalHandler.hpp"
#include "utils/EventUtils.hpp"
#include <boost/filesystem.hpp>


namespace precice {
namespace utils {

void terminationSignalHandler(int signal)
{
  // Print the events statistics
  precice::utils::EventRegistry::instance().signal_handler(signal);

  // Remove connection info files. This is just a guess and will only
  // work if the address directory is ".".
  boost::filesystem::remove_all("precice-run");
}

}
}

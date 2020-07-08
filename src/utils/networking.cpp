#include "utils/networking.hpp"

namespace precice {
namespace utils {
namespace networking {

std::string loopbackInterfaceName()
{
#if defined(__linux__)
  return "lo";
#elif defined(__APPLE__)
  return "lo0";
#elif defined(__WIN32)
  // Not required as we directly use the 127.0.0.1 under Windows
  return "";
#else
#warning "Your target architecture does not define a loopback interface. Please consider reporting this to the preCICE developers."
  return "";
#endif
}

} // namespace networking
} // namespace utils
} // namespace precice

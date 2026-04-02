#include "utils/networking.hpp"

namespace precice::utils::networking {

std::string loopbackInterfaceName()
{
#if defined(__linux__)
  return "lo";
#elif defined(__APPLE__) || defined(BSD)
  return "lo0";
#elif defined(_WIN32)
  // Not required as we directly use the 127.0.0.1 under Windows
  return "";
#else
#error "There is no loopback device defined for your OS. Please open an issue with a suggestion at https://github.com/precice/precice/issues"
#endif
}

} // namespace precice::utils::networking

#pragma once

#include <string>

namespace precice {
namespace utils {
namespace networking {

/// Returns the name of the canonical loopback interface on this system
std::string loopbackInterfaceName();

} // namespace networking

} // namespace utils
} // namespace precice

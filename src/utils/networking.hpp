#pragma once

#include <string>

namespace precice::utils::networking {

/// Returns the name of the canonical loopback interface on this system
std::string loopbackInterfaceName();

} // namespace precice::utils::networking

#pragma once

#include <string>
#include <string_view>

namespace precice::utils {

/// creates a portable hash of the given input
std::string preciceHash(std::string_view s);

} // namespace precice::utils

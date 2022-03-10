#pragma once

#include <stdexcept>

namespace precice {

class Error : public std::runtime_error {
public:
  Error(const std::string &what_arg)
      : std::runtime_error(what_arg){};
};

} // namespace precice

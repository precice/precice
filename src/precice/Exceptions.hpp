#pragma once

#include <stdexcept>

#include "precice/export.h"

namespace precice {

class PRECICE_API Error : public std::runtime_error {
public:
  Error(const std::string &what_arg)
      : std::runtime_error(what_arg){};
};

} // namespace precice

#pragma once

#include <vector>
#include "com/SharedPointer.hpp"

namespace precice::com {
class Request {

public:
  static void wait(std::vector<PtrRequest> &requests);

  virtual ~Request();

  virtual bool test() = 0;

  virtual void wait() = 0;
};
} // namespace precice::com

#pragma once

#include "com/SharedPointer.hpp"
#include <vector>

namespace precice {
namespace com {
class Request {

public:
  static void wait(std::vector<PtrRequest>& requests);

  virtual ~Request();

  virtual bool test() = 0;

  virtual void wait() = 0;
};

}} // namespace precice, com

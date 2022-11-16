#include <memory>

#include "com/Request.hpp"

namespace precice::com {

void Request::wait(std::vector<PtrRequest> &requests)
{
  for (const auto &request : requests) {
    request->wait();
  }
}

Request::~Request() = default;
} // namespace precice::com

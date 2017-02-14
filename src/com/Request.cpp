#include "Request.hpp"

namespace precice {
namespace com {

void Request::wait(std::vector<PtrRequest>& requests) {
  for (auto request : requests) {
    request->wait();
  }
}

Request::~Request() {
}
}
} // namespace precice, com

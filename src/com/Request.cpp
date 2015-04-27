// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "Request.hpp"

namespace precice {
namespace com {
void
Request::wait(std::vector<SharedPointer>& requests) {
  for (auto request : requests) {
    request->wait();
  }
}

Request::~Request() {
}
}
} // namespace precice, com

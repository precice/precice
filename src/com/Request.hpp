// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_REQUEST_HPP_
#define PRECICE_COM_REQUEST_HPP_

#include <memory>
#include <vector>

namespace precice {
namespace com {
class Request {
public:
  using SharedPointer = std::shared_ptr<Request>;

public:
  static void wait(std::vector<SharedPointer>& requests);

  virtual ~Request();

  virtual bool test() = 0;

  virtual void wait() = 0;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_REQUEST_HPP_ */

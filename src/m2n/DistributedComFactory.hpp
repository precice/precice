#pragma once

#include "DistributedCommunication.hpp"

#include <memory>

namespace precice::m2n {
class DistributedComFactory {

public:
  using SharedPointer = std::shared_ptr<DistributedComFactory>;

  virtual ~DistributedComFactory() = default;

  virtual DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) = 0;
};
} // namespace precice::m2n

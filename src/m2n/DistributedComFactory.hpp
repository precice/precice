#pragma once

#include "DistributedCommunication.hpp"

#include <memory>

namespace precice {
namespace m2n {
class DistributedComFactory {
public:
  using SharedPointer = std::shared_ptr<DistributedComFactory>;

public:
  /**
   * @brief Destructor.
   */
  virtual ~DistributedComFactory(){};

  virtual DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) = 0;
};
}
} // namespace precice, m2n

#pragma once

#include "DistributedComFactory.hpp"

namespace precice {
namespace m2n {
class GatherScatterComFactory : public DistributedComFactory {
public:
  /**
   * @brief Constructor.
   */
  GatherScatterComFactory(com::Communication::SharedPointer masterCom);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh);

private:
  // @brief communication between the master processes
  com::Communication::SharedPointer _masterCom;
};
}
} // namespace precice, m2n

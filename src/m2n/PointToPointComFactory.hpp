#pragma once

#include "DistributedComFactory.hpp"

#include "com/CommunicationFactory.hpp"

namespace precice {
namespace m2n {
class PointToPointComFactory : public DistributedComFactory {
public:
  /**
   * @brief Constructor.
   */
  PointToPointComFactory(com::CommunicationFactory::SharedPointer comFactory);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh);

private:
  // @brief communication factory for 1:M communications
  com::CommunicationFactory::SharedPointer _comFactory;
};
}
} // namespace precice, m2n

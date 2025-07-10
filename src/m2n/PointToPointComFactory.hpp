#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::m2n {

class PointToPointComFactory : public DistributedComFactory {

public:
  explicit PointToPointComFactory(com::PtrCommunicationFactory comFactory);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) override;

private:
  /// communication factory for 1:M communications
  com::PtrCommunicationFactory _comFactory;
};

} // namespace precice::m2n

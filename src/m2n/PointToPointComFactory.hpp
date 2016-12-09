#pragma once

#include "DistributedComFactory.hpp"

#include "com/CommunicationFactory.hpp"

namespace precice {
namespace m2n {

class PointToPointComFactory : public DistributedComFactory {

public:
  explicit PointToPointComFactory(com::CommunicationFactory::SharedPointer comFactory);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh);

private:
  /// communication factory for 1:M communications
  com::CommunicationFactory::SharedPointer _comFactory;

};


}} // namespace precice, m2n

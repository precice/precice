#include "PointToPointComFactory.hpp"
#include "PointToPointCommunication.hpp"

#include "com/SharedPointer.hpp"

namespace precice {
namespace m2n {

PointToPointComFactory::PointToPointComFactory(com::PtrCommunicationFactory comFactory)
    : _comFactory(comFactory) {}

DistributedCommunication::SharedPointer
PointToPointComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(new PointToPointCommunication(_comFactory, mesh));
}

} // namespace m2n
} // namespace precice

#include "PointToPointCommunication.hpp"

#include "PointToPointComFactory.hpp"

namespace precice {
namespace m2n {

PointToPointComFactory::PointToPointComFactory(com::CommunicationFactory::SharedPointer comFactory)
  : _comFactory(comFactory) {}

DistributedCommunication::SharedPointer
PointToPointComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(new PointToPointCommunication(_comFactory, mesh));
}

}} // namespace precice, m2n

#include "precice/impl/DataContext.hpp"
#include <memory>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  _toData   = data;
  _fromData = data;
  _mesh     = mesh;
}

std::string DataContext::getDataName() const
{
  return _fromData->getName();
}

std::string DataContext::getMeshName() const
{
  return _mesh->getName();
}

int DataContext::getMeshID() const
{
  return _mesh->getID();
}

const mesh::PtrData DataContext::fromData() const
{
  return _fromData;
}

void DataContext::setFromData(mesh::PtrData data)
{
  _fromData = data;
}

const mesh::PtrData DataContext::toData() const
{
  return _toData;
}

void DataContext::setToData(mesh::PtrData data)
{
  _toData = data;
}

} // namespace impl
} // namespace precice

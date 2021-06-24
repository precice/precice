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

std::string DataContext::getFromDataName() const
{
  return _fromData->getName();
}

int DataContext::getFromDataID() const
{
  return _fromData->getID();
}

std::string DataContext::getToDataName() const
{
  return _toData->getName();
}

int DataContext::getToDataID() const
{
  return _toData->getID();
}

std::string DataContext::getMeshName() const
{
  return _mesh->getName();
}

int DataContext::getMeshID() const
{
  return _mesh->getID();
}

mesh::PtrData DataContext::fromData()
{
  return _fromData;
}

void DataContext::setFromData(mesh::PtrData data)
{
  _fromData = data;
}

mesh::PtrData DataContext::toData()
{
  return _toData;
}

void DataContext::setToData(mesh::PtrData data)
{
  _toData = data;
}

} // namespace impl
} // namespace precice

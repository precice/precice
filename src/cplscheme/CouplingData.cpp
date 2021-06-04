#include "cplscheme/CouplingData.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

CouplingData::CouplingData(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
    : requiresInitialization(requiresInitialization),
      _data(data),
      _mesh(mesh)
{
  PRECICE_ASSERT(_data != nullptr);
  PRECICE_ASSERT(_mesh != nullptr);
  PRECICE_ASSERT(_mesh.use_count() > 0);
}

int CouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
}

Eigen::VectorXd &CouplingData::values()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

const Eigen::VectorXd &CouplingData::values() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

void CouplingData::storeIteration()
{
  _previousIteration = this->values();
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousIteration;
}

int CouplingData::getMeshID()
{
  return _mesh->getID();
}

int CouplingData::getDataID()
{
  return _data->getID();
}

std::string CouplingData::getDataName()
{
  return _data->getName();
}

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

} // namespace cplscheme
} // namespace precice
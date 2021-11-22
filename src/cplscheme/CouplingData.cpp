#include "cplscheme/CouplingData.hpp"

#include <utility>

#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "time/Waveform.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

CouplingData::CouplingData(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
    : requiresInitialization(requiresInitialization),
      _data(std::move(data)),
      _mesh(std::move(mesh))
{
  PRECICE_ASSERT(_data != nullptr);
  _previousIteration = Eigen::VectorXd::Zero(_data->values().size());
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

void CouplingData::initializeStorage()
{
  _data->waveform()->initialize(values().size());
  storeIteration();
}

void CouplingData::moveToNextWindow()
{
  _data->waveform()->moveToNextWindow();
  values() = _data->waveform()->getInitialGuess();
}

void CouplingData::storeDataInWaveform()
{
  _data->waveform()->store(values());
}

} // namespace cplscheme
} // namespace precice
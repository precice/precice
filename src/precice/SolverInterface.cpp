#include "precice/SolverInterface.hpp"
#include "cplscheme/Constants.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/versions.hpp"

namespace precice {

SolverInterface::SolverInterface(
    const std::string &participantName,
    const std::string &configurationFileName,
    int                solverProcessIndex,
    int                solverProcessSize)
    : _impl(new impl::SolverInterfaceImpl(participantName, configurationFileName, solverProcessIndex, solverProcessSize))
{
}

SolverInterface::SolverInterface(
    const std::string &participantName,
    const std::string &configurationFileName,
    int                solverProcessIndex,
    int                solverProcessSize,
    void *             communicator)
    : _impl(new impl::SolverInterfaceImpl(participantName, configurationFileName, solverProcessIndex, solverProcessSize, communicator))
{
}

SolverInterface::~SolverInterface() = default;

double SolverInterface::initialize()
{
  return _impl->initialize();
}

void SolverInterface::initializeData()
{
  _impl->initializeData();
}

double SolverInterface::advance(
    double computedTimestepLength)
{
  return _impl->advance(computedTimestepLength);
}

void SolverInterface::finalize()
{
  return _impl->finalize();
}

int SolverInterface::getDimensions() const
{
  return _impl->getDimensions();
}

bool SolverInterface::isCouplingOngoing() const
{
  return _impl->isCouplingOngoing();
}

bool SolverInterface::isReadDataAvailable() const
{
  return _impl->isReadDataAvailable();
}

bool SolverInterface::isWriteDataRequired(
    double computedTimestepLength) const
{
  return _impl->isWriteDataRequired(computedTimestepLength);
}

bool SolverInterface::isTimeWindowComplete() const
{
  return _impl->isTimeWindowComplete();
}

bool SolverInterface::isActionRequired(
    const std::string &action) const
{
  return _impl->isActionRequired(action);
}

void SolverInterface::markActionFulfilled(
    const std::string &action)
{
  _impl->markActionFulfilled(action);
}

bool SolverInterface::hasMesh(
    const std::string &meshName) const
{
  return _impl->hasMesh(meshName);
}

int SolverInterface::getMeshID(
    const std::string &meshName) const
{
  return _impl->getMeshID(meshName);
}

std::set<int> SolverInterface::getMeshIDs() const
{
  return _impl->getMeshIDs();
}

bool SolverInterface::hasData(
    const std::string &dataName, int meshID) const
{
  return _impl->hasData(dataName, meshID);
}

int SolverInterface::getDataID(
    const std::string &dataName, int meshID) const
{
  return _impl->getDataID(dataName, meshID);
}

bool SolverInterface::hasToEvaluateSurrogateModel() const
{
  return _impl->hasToEvaluateSurrogateModel();
}

bool SolverInterface::hasToEvaluateFineModel() const
{
  return _impl->hasToEvaluateFineModel();
}

//void SolverInterface:: resetMesh
//(
//  int meshID )
//{
//  _impl->resetMesh(meshID);
//}

int SolverInterface::setMeshVertex(
    int           meshID,
    const double *position)
{
  return _impl->setMeshVertex(meshID, position);
}

int SolverInterface::getMeshVertexSize(
    int meshID) const
{
  return _impl->getMeshVertexSize(meshID);
}

void SolverInterface::setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  _impl->setMeshVertices(meshID, size, positions, ids);
}

void SolverInterface::getMeshVertices(
    int        meshID,
    int        size,
    const int *ids,
    double *   positions) const
{
  _impl->getMeshVertices(meshID, size, ids, positions);
}

void SolverInterface::getMeshVertexIDsFromPositions(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids) const
{
  _impl->getMeshVertexIDsFromPositions(meshID, size, positions, ids);
}

int SolverInterface::setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
{
  return _impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void SolverInterface::setMeshTriangle(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID)
{
  _impl->setMeshTriangle(meshID, firstEdgeID, secondEdgeID, thirdEdgeID);
}

void SolverInterface::setMeshTriangleWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
{
  _impl->setMeshTriangleWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void SolverInterface::setMeshQuad(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID)
{
  _impl->setMeshQuad(meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID);
}

void SolverInterface::setMeshQuadWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  _impl->setMeshQuadWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID,
                              fourthVertexID);
}

void SolverInterface::mapReadDataTo(
    int toMeshID)
{
  _impl->mapReadDataTo(toMeshID);
}

void SolverInterface::mapWriteDataFrom(
    int fromMeshID)
{
  _impl->mapWriteDataFrom(fromMeshID);
}

void SolverInterface::writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  _impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *value)
{
  _impl->writeVectorData(dataID, valueIndex, value);
}

void SolverInterface::writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  _impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::writeScalarData(
    int    dataID,
    int    valueIndex,
    double value)
{
  _impl->writeScalarData(dataID, valueIndex, value);
}

void SolverInterface::readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  _impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::readVectorData(
    int     dataID,
    int     valueIndex,
    double *value) const
{
  return _impl->readVectorData(dataID, valueIndex, value);
}

void SolverInterface::readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  _impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::readScalarData(
    int     dataID,
    int     valueIndex,
    double &value) const
{
  return _impl->readScalarData(dataID, valueIndex, value);
}

std::string getVersionInformation()
{
  return {precice::versionInformation};
}

namespace constants {

const std::string &actionWriteInitialData()
{
  return cplscheme::constants::actionWriteInitialData();
}

const std::string &actionWriteIterationCheckpoint()
{
  return cplscheme::constants::actionWriteIterationCheckpoint();
}

const std::string &actionReadIterationCheckpoint()
{
  return cplscheme::constants::actionReadIterationCheckpoint();
}

} // namespace constants

} // namespace precice

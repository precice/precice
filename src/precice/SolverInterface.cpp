#include "precice/SolverInterface.hpp"
#include "cplscheme/Constants.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/versions.hpp"
#include "precice/types.hpp"

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

MeshID SolverInterface::getMeshID(
    const std::string &meshName) const
{
  return _impl->getMeshID(meshName);
}

std::set<int> SolverInterface::getMeshIDs() const
{
  return _impl->getMeshIDs();
}

bool SolverInterface::isMeshConnectivityRequired(MeshID meshID) const
{
  return _impl->isMeshConnectivityRequired(meshID);
}

bool SolverInterface::isGradientDataRequired(DataID dataID) const
{
  return _impl->isGradientDataRequired(dataID);
}

bool SolverInterface::hasData(
    const std::string &dataName, MeshID meshID) const
{
  return _impl->hasData(dataName, meshID);
}

DataID SolverInterface::getDataID(
    const std::string &dataName, MeshID meshID) const
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
//  MeshID meshID )
//{
//  _impl->resetMesh(meshID);
//}

VertexID SolverInterface::setMeshVertex(
    MeshID        meshID,
    const double *position)
{
  return _impl->setMeshVertex(meshID, position);
}

Size SolverInterface::getMeshVertexSize(
    MeshID meshID) const
{
  return _impl->getMeshVertexSize(meshID);
}

void SolverInterface::setMeshVertices(
    MeshID        meshID,
    Size          size,
    const double *positions,
    VertexID *    ids)
{
  _impl->setMeshVertices(meshID, size, positions, ids);
}

void SolverInterface::getMeshVertices(
    MeshID          meshID,
    Size            size,
    const VertexID *ids,
    double *        positions) const
{
  _impl->getMeshVertices(meshID, size, ids, positions);
}

void SolverInterface::getMeshVertexIDsFromPositions(
    MeshID        meshID,
    Size          size,
    const double *positions,
    VertexID *    ids) const
{
  _impl->getMeshVertexIDsFromPositions(meshID, size, positions, ids);
}

EdgeID SolverInterface::setMeshEdge(
    MeshID   meshID,
    VertexID firstVertexID,
    VertexID secondVertexID)
{
  return _impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void SolverInterface::setMeshTriangle(
    MeshID meshID,
    EdgeID firstEdgeID,
    EdgeID secondEdgeID,
    EdgeID thirdEdgeID)
{
  _impl->setMeshTriangle(meshID, firstEdgeID, secondEdgeID, thirdEdgeID);
}

void SolverInterface::setMeshTriangleWithEdges(
    MeshID   meshID,
    VertexID firstVertexID,
    VertexID secondVertexID,
    VertexID thirdVertexID)
{
  _impl->setMeshTriangleWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void SolverInterface::setMeshQuad(
    MeshID meshID,
    EdgeID firstEdgeID,
    EdgeID secondEdgeID,
    EdgeID thirdEdgeID,
    EdgeID fourthEdgeID)
{
  _impl->setMeshQuad(meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID);
}

void SolverInterface::setMeshQuadWithEdges(
    MeshID   meshID,
    VertexID firstVertexID,
    VertexID secondVertexID,
    VertexID thirdVertexID,
    VertexID fourthVertexID)
{
  _impl->setMeshQuadWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID,
                              fourthVertexID);
}

void SolverInterface::mapReadDataTo(
    MeshID toMeshID)
{
  _impl->mapReadDataTo(toMeshID);
}

void SolverInterface::mapWriteDataFrom(
    MeshID fromMeshID)
{
  _impl->mapWriteDataFrom(fromMeshID);
}

void SolverInterface::writeBlockVectorData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    const double *  values)
{
  _impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::writeBlockVectorGradientData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    const double *  gradientValues,
    bool            rowsFirst)
{
  _impl->writeBlockVectorGradientData(dataID, size, valueIndices, gradientValues, rowsFirst);
}

void SolverInterface::writeVectorData(
    DataID        dataID,
    VertexID      valueIndex,
    const double *value)
{
  _impl->writeVectorData(dataID, valueIndex, value);
}

void SolverInterface::writeVectorGradientData(
    DataID        dataID,
    VertexID      valueIndex,
    const double *gradientValues,
    bool          rowsFirst)
{
  _impl->writeVectorGradientData(dataID, valueIndex, gradientValues, rowsFirst);
}

void SolverInterface::writeBlockScalarData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    const double *  values)
{
  _impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::writeBlockScalarGradientData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    const double *  gradientValues)
{
  _impl->writeBlockScalarGradientData(dataID, size, valueIndices, gradientValues);
}

void SolverInterface::writeScalarData(
    DataID   dataID,
    VertexID valueIndex,
    double   value)
{
  _impl->writeScalarData(dataID, valueIndex, value);
}

void SolverInterface::writeScalarGradientData(
    DataID        dataID,
    VertexID      valueIndex,
    const double *gradientValues)
{
  _impl->writeScalarGradientData(dataID, valueIndex, gradientValues);
}

void SolverInterface::readBlockVectorData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    double *        values) const
{
  _impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::readBlockVectorData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    double          relativeReadTime,
    double *        values) const
{
  _impl->readBlockVectorData(dataID, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readVectorData(
    DataID   dataID,
    VertexID valueIndex,
    double * value) const
{
  _impl->readVectorData(dataID, valueIndex, value);
}

void SolverInterface::readVectorData(
    DataID   dataID,
    VertexID valueIndex,
    double   relativeReadTime,
    double * value) const
{
  // @todo: needs testing!
  _impl->readVectorData(dataID, valueIndex, relativeReadTime, value);
}

void SolverInterface::readBlockScalarData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    double *        values) const
{
  _impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::readBlockScalarData(
    DataID          dataID,
    Size            size,
    const VertexID *valueIndices,
    double          relativeReadTime,
    double *        values) const
{
  _impl->readBlockScalarData(dataID, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readScalarData(
    DataID   dataID,
    VertexID valueIndex,
    double & value) const
{
  _impl->readScalarData(dataID, valueIndex, value);
}

void SolverInterface::readScalarData(
    DataID   dataID,
    VertexID valueIndex,
    double   relativeReadTime,
    double & value) const
{
  _impl->readScalarData(dataID, valueIndex, relativeReadTime, value);
}

void SolverInterface::setMeshAccessRegion(const MeshID  meshID,
                                          const double *boundingBox) const
{
  _impl->setMeshAccessRegion(meshID, boundingBox);
}

void SolverInterface::getMeshVerticesAndIDs(const MeshID meshID,
                                            const Size   size,
                                            VertexID *   ids,
                                            double *     coordinates) const
{
  _impl->getMeshVerticesAndIDs(meshID, size, ids, coordinates);
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

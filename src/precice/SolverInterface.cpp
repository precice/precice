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

bool SolverInterface::isTimeWindowComplete() const
{
  return _impl->isTimeWindowComplete();
}

bool SolverInterface::requiresInitialData()
{
  return _impl->requiresInitialData();
}

bool SolverInterface::requiresReadingCheckpoint()
{
  return _impl->requiresReadingCheckpoint();
}

bool SolverInterface::requiresWritingCheckpoint()
{
  return _impl->requiresWritingCheckpoint();
}

bool SolverInterface::hasMesh(
    const std::string &meshName) const
{
  return _impl->hasMesh(meshName);
}

bool SolverInterface::requiresMeshConnectivityFor(std::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(meshName);
}

bool SolverInterface::requiresGradientDataFor(std::string_view meshName,
                                              std::string_view dataName) const
{
  return _impl->requiresGradientDataFor(meshName, dataName);
}

bool SolverInterface::hasData(std::string_view meshName, std::string_view dataName) const
{
  return _impl->hasData(meshName, dataName);
}

// void SolverInterface:: resetMesh
//(
//   std::string_view meshName )
//{
//   _impl->resetMesh(meshNameID);
// }

int SolverInterface::setMeshVertex(
    std::string_view meshName,
    const double *   position)
{
  return _impl->setMeshVertex(meshName, position);
}

int SolverInterface::getMeshVertexSize(
    std::string_view meshName) const
{
  return _impl->getMeshVertexSize(meshName);
}

void SolverInterface::setMeshVertices(
    std::string_view meshName,
    int              size,
    const double *   positions,
    int *            ids)
{
  _impl->setMeshVertices(meshName, size, positions, ids);
}

void SolverInterface::setMeshEdge(
    std::string_view meshName,
    int              firstVertexID,
    int              secondVertexID)
{
  _impl->setMeshEdge(meshName, firstVertexID, secondVertexID);
}

void SolverInterface::setMeshEdges(
    std::string_view meshName,
    int              size,
    const int *      vertices)
{
  _impl->setMeshEdges(meshName, size, vertices);
}

void SolverInterface::setMeshTriangle(
    std::string_view meshName,
    int              firstVertexID,
    int              secondVertexID,
    int              thirdVertexID)
{
  _impl->setMeshTriangle(meshName, firstVertexID, secondVertexID, thirdVertexID);
}

void SolverInterface::setMeshTriangles(
    std::string_view meshName,
    int              size,
    const int *      vertices)
{
  _impl->setMeshTriangles(meshName, size, vertices);
}

void SolverInterface::setMeshQuad(
    std::string_view meshName,
    int              firstVertexID,
    int              secondVertexID,
    int              thirdVertexID,
    int              fourthVertexID)
{
  _impl->setMeshQuad(meshName, firstVertexID, secondVertexID, thirdVertexID,
                     fourthVertexID);
}

void SolverInterface::setMeshQuads(
    std::string_view meshName,
    int              size,
    const int *      vertices)
{
  _impl->setMeshQuads(meshName, size, vertices);
}

void SolverInterface::setMeshTetrahedron(
    std::string_view meshName,
    int              firstVertexID,
    int              secondVertexID,
    int              thirdVertexID,
    int              fourthVertexID)
{
  _impl->setMeshTetrahedron(meshName, firstVertexID, secondVertexID, thirdVertexID,
                            fourthVertexID);
}

void SolverInterface::setMeshTetrahedra(
    std::string_view meshName,
    int              size,
    const int *      vertices)
{
  _impl->setMeshTetrahedra(meshName, size, vertices);
}

void SolverInterface::writeBlockVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    const double *   values)
{
  _impl->writeBlockVectorData(meshName, dataName, size, valueIndices, values);
}

void SolverInterface::writeBlockVectorGradientData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    const double *   gradientValues)
{
  _impl->writeBlockVectorGradientData(meshName, dataName, size, valueIndices, gradientValues);
}

void SolverInterface::writeVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    const double *   value)
{
  _impl->writeVectorData(meshName, dataName, valueIndex, value);
}

void SolverInterface::writeVectorGradientData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    const double *   gradientValues)
{
  _impl->writeVectorGradientData(meshName, dataName, valueIndex, gradientValues);
}

void SolverInterface::writeBlockScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    const double *   values)
{
  _impl->writeBlockScalarData(meshName, dataName, size, valueIndices, values);
}

void SolverInterface::writeBlockScalarGradientData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    const double *   gradientValues)
{
  _impl->writeBlockScalarGradientData(meshName, dataName, size, valueIndices, gradientValues);
}

void SolverInterface::writeScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    double           value)
{
  _impl->writeScalarData(meshName, dataName, valueIndex, value);
}

void SolverInterface::writeScalarGradientData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    const double *   gradientValues)
{
  _impl->writeScalarGradientData(meshName, dataName, valueIndex, gradientValues);
}

void SolverInterface::readBlockVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    double *         values) const
{
  _impl->readBlockVectorData(meshName, dataName, size, valueIndices, values);
}

void SolverInterface::readBlockVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    double           relativeReadTime,
    double *         values) const
{
  _impl->readBlockVectorData(meshName, dataName, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    double *         value) const
{
  _impl->readVectorData(meshName, dataName, valueIndex, value);
}

void SolverInterface::readVectorData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    double           relativeReadTime,
    double *         value) const
{
  // @todo: needs testing!
  _impl->readVectorData(meshName, dataName, valueIndex, relativeReadTime, value);
}

void SolverInterface::readBlockScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    double *         values) const
{
  _impl->readBlockScalarData(meshName, dataName, size, valueIndices, values);
}

void SolverInterface::readBlockScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              size,
    const int *      valueIndices,
    double           relativeReadTime,
    double *         values) const
{
  _impl->readBlockScalarData(meshName, dataName, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    double &         value) const
{
  _impl->readScalarData(meshName, dataName, valueIndex, value);
}

void SolverInterface::readScalarData(
    std::string_view meshName,
    std::string_view dataName,
    int              valueIndex,
    double           relativeReadTime,
    double &         value) const
{
  _impl->readScalarData(meshName, dataName, valueIndex, relativeReadTime, value);
}

void SolverInterface::setMeshAccessRegion(std::string_view meshName,
                                          const double *   boundingBox) const
{
  _impl->setMeshAccessRegion(meshName, boundingBox);
}

void SolverInterface::getMeshVerticesAndIDs(std::string_view meshName,
                                            const int        size,
                                            int *            ids,
                                            double *         coordinates) const
{
  _impl->getMeshVerticesAndIDs(meshName, size, ids, coordinates);
}

} // namespace precice

#include "precice/SolverInterface.hpp"
#include <string_view>
#include "cplscheme/Constants.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/versions.hpp"

namespace precice {

namespace {
std::string_view toSV(precice::string_view sv)
{
  return {sv.data(), sv.size()};
}
} // namespace

SolverInterface::SolverInterface(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize)
    : _impl(new impl::SolverInterfaceImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize))
{
}

SolverInterface::SolverInterface(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize,
    void *                 communicator)
    : _impl(new impl::SolverInterfaceImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, communicator))
{
}

SolverInterface::~SolverInterface() = default;

void SolverInterface::initialize()
{
  _impl->initialize();
}

void SolverInterface::advance(
    double computedTimeStepSize)
{
  _impl->advance(computedTimeStepSize);
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

double SolverInterface::getMaxTimeStepSize() const
{
  return _impl->getMaxTimeStepSize();
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

bool SolverInterface::hasMesh(::precice::string_view meshName) const
{
  return _impl->hasMesh(toSV(toSV(meshName)));
}

bool SolverInterface::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(toSV(toSV(meshName)));
}

bool SolverInterface::requiresGradientDataFor(::precice::string_view meshName,
                                              ::precice::string_view dataName) const
{
  return _impl->requiresGradientDataFor(toSV(meshName), toSV(dataName));
}

bool SolverInterface::hasData(::precice::string_view meshName, ::precice::string_view dataName) const
{
  return _impl->hasData(toSV(meshName), toSV(dataName));
}

// void SolverInterface:: resetMesh
//(
//   ::precice::string_view meshName )
//{
//   _impl->resetMesh(toSV(meshName)ID);
// }

int SolverInterface::setMeshVertex(
    ::precice::string_view meshName,
    const double *         position)
{
  return _impl->setMeshVertex(toSV(meshName), position);
}

int SolverInterface::getMeshVertexSize(
    ::precice::string_view meshName) const
{
  return _impl->getMeshVertexSize(toSV(meshName));
}

void SolverInterface::setMeshVertices(
    ::precice::string_view meshName,
    int                    size,
    const double *         positions,
    int *                  ids)
{
  _impl->setMeshVertices(toSV(meshName), size, positions, ids);
}

void SolverInterface::setMeshEdge(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID)
{
  _impl->setMeshEdge(toSV(meshName), firstVertexID, secondVertexID);
}

void SolverInterface::setMeshEdges(
    ::precice::string_view meshName,
    int                    size,
    const int *            vertices)
{
  _impl->setMeshEdges(toSV(meshName), size, vertices);
}

void SolverInterface::setMeshTriangle(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID)
{
  _impl->setMeshTriangle(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID);
}

void SolverInterface::setMeshTriangles(
    ::precice::string_view meshName,
    int                    size,
    const int *            vertices)
{
  _impl->setMeshTriangles(toSV(meshName), size, vertices);
}

void SolverInterface::setMeshQuad(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID,
    int                    fourthVertexID)
{
  _impl->setMeshQuad(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID,
                     fourthVertexID);
}

void SolverInterface::setMeshQuads(
    ::precice::string_view meshName,
    int                    size,
    const int *            vertices)
{
  _impl->setMeshQuads(toSV(meshName), size, vertices);
}

void SolverInterface::setMeshTetrahedron(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID,
    int                    fourthVertexID)
{
  _impl->setMeshTetrahedron(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID,
                            fourthVertexID);
}

void SolverInterface::setMeshTetrahedra(
    ::precice::string_view meshName,
    int                    size,
    const int *            vertices)
{
  _impl->setMeshTetrahedra(toSV(meshName), size, vertices);
}

void SolverInterface::writeBlockVectorData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    const double *         values)
{
  _impl->writeBlockVectorData(toSV(meshName), toSV(dataName), size, valueIndices, values);
}

void SolverInterface::writeBlockVectorGradientData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    const double *         gradientValues)
{
  _impl->writeBlockVectorGradientData(toSV(meshName), toSV(dataName), size, valueIndices, gradientValues);
}

void SolverInterface::writeVectorData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    const double *         value)
{
  _impl->writeVectorData(toSV(meshName), toSV(dataName), valueIndex, value);
}

void SolverInterface::writeVectorGradientData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    const double *         gradientValues)
{
  _impl->writeVectorGradientData(toSV(meshName), toSV(dataName), valueIndex, gradientValues);
}

void SolverInterface::writeBlockScalarData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    const double *         values)
{
  _impl->writeBlockScalarData(toSV(meshName), toSV(dataName), size, valueIndices, values);
}

void SolverInterface::writeBlockScalarGradientData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    const double *         gradientValues)
{
  _impl->writeBlockScalarGradientData(toSV(meshName), toSV(dataName), size, valueIndices, gradientValues);
}

void SolverInterface::writeScalarData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    double                 value)
{
  _impl->writeScalarData(toSV(meshName), toSV(dataName), valueIndex, value);
}

void SolverInterface::writeScalarGradientData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    const double *         gradientValues)
{
  _impl->writeScalarGradientData(toSV(meshName), toSV(dataName), valueIndex, gradientValues);
}

void SolverInterface::readBlockVectorData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    double                 relativeReadTime,
    double *               values) const
{
  _impl->readBlockVectorData(toSV(meshName), toSV(dataName), size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readVectorData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    double                 relativeReadTime,
    double *               value) const
{
  // @todo: needs testing!
  _impl->readVectorData(toSV(meshName), toSV(dataName), valueIndex, relativeReadTime, value);
}

void SolverInterface::readBlockScalarData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    size,
    const int *            valueIndices,
    double                 relativeReadTime,
    double *               values) const
{
  _impl->readBlockScalarData(toSV(meshName), toSV(dataName), size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readScalarData(
    ::precice::string_view meshName,
    ::precice::string_view dataName,
    int                    valueIndex,
    double                 relativeReadTime,
    double &               value) const
{
  _impl->readScalarData(toSV(meshName), toSV(dataName), valueIndex, relativeReadTime, value);
}

void SolverInterface::setMeshAccessRegion(::precice::string_view meshName,
                                          const double *         boundingBox) const
{
  _impl->setMeshAccessRegion(toSV(meshName), boundingBox);
}

void SolverInterface::getMeshVerticesAndIDs(::precice::string_view meshName,
                                            const int              size,
                                            int *                  ids,
                                            double *               coordinates) const
{
  _impl->getMeshVerticesAndIDs(toSV(meshName), size, ids, coordinates);
}

} // namespace precice

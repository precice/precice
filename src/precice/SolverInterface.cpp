#include "precice/SolverInterface.hpp"
#include <string_view>
#include "cplscheme/Constants.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/versions.hpp"

namespace precice {

namespace {
std::string_view toSV(precice::string_view sv)
{
  // The given string_view may contain null chars.
  // We trim to the first null char here.
  std::string_view s{sv.data(), sv.size()};
  auto             trim_pos = s.find('\0');
  if (trim_pos != s.npos) {
    s.remove_suffix(s.size() - trim_pos);
  }
  return s;
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

int SolverInterface::getMeshDimensions(::precice::string_view meshName) const
{
  return _impl->getMeshDimensions(toSV(meshName));
}

int SolverInterface::getDataDimensions(::precice::string_view meshName, ::precice::string_view dataName) const
{
  return _impl->getDataDimensions(toSV(meshName), toSV(dataName));
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
  return _impl->hasMesh(toSV(meshName));
}

bool SolverInterface::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(toSV(meshName));
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

void SolverInterface::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   values)
{
  _impl->writeData(toSV(meshName), toSV(dataName), vertices, values);
}

void SolverInterface::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  _impl->readData(toSV(meshName), toSV(dataName), vertices, relativeReadTime, values);
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

void SolverInterface::writeGradientData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   gradients)
{
  _impl->writeGradientData(toSV(meshName), toSV(dataName), vertices, gradients);
}

} // namespace precice

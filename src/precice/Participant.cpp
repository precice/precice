#include <optional>
#include <string_view>

#include "cplscheme/Constants.hpp"
#include "precice/Participant.hpp"
#include "precice/impl/ParticipantImpl.hpp"
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

Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize)
    : _impl(new impl::ParticipantImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, std::nullopt))
{
}

Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize,
    void *                 communicator)
    : _impl(new impl::ParticipantImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, {communicator}))
{
}

Participant::~Participant() = default;

void Participant::initialize()
{
  _impl->initialize();
}

void Participant::advance(
    double computedTimeStepSize)
{
  _impl->advance(computedTimeStepSize);
}

void Participant::finalize()
{
  return _impl->finalize();
}

int Participant::getMeshDimensions(::precice::string_view meshName) const
{
  return _impl->getMeshDimensions(toSV(meshName));
}

int Participant::getDataDimensions(::precice::string_view meshName, ::precice::string_view dataName) const
{
  return _impl->getDataDimensions(toSV(meshName), toSV(dataName));
}

bool Participant::isCouplingOngoing() const
{
  return _impl->isCouplingOngoing();
}

bool Participant::isTimeWindowComplete() const
{
  return _impl->isTimeWindowComplete();
}

double Participant::getMaxTimeStepSize() const
{
  return _impl->getMaxTimeStepSize();
}

bool Participant::requiresInitialData()
{
  return _impl->requiresInitialData();
}

bool Participant::requiresReadingCheckpoint()
{
  return _impl->requiresReadingCheckpoint();
}

bool Participant::requiresWritingCheckpoint()
{
  return _impl->requiresWritingCheckpoint();
}

bool Participant::hasMesh(::precice::string_view meshName) const
{
  return _impl->hasMesh(toSV(meshName));
}

bool Participant::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(toSV(meshName));
}

bool Participant::requiresGradientDataFor(::precice::string_view meshName,
                                          ::precice::string_view dataName) const
{
  return _impl->requiresGradientDataFor(toSV(meshName), toSV(dataName));
}

bool Participant::hasData(::precice::string_view meshName, ::precice::string_view dataName) const
{
  return _impl->hasData(toSV(meshName), toSV(dataName));
}

// void SolverInterface:: resetMesh
//(
//   ::precice::string_view meshName )
//{
//   _impl->resetMesh(toSV(meshName)ID);
// }

int Participant::setMeshVertex(
    ::precice::string_view        meshName,
    ::precice::span<const double> position)
{
  return _impl->setMeshVertex(toSV(meshName), position);
}

int Participant::getMeshVertexSize(
    ::precice::string_view meshName) const
{
  return _impl->getMeshVertexSize(toSV(meshName));
}

void Participant::setMeshVertices(
    ::precice::string_view        meshName,
    ::precice::span<const double> positions,
    ::precice::span<VertexID>     ids)
{
  _impl->setMeshVertices(toSV(meshName), positions, ids);
}

void Participant::setMeshEdge(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID)
{
  _impl->setMeshEdge(toSV(meshName), firstVertexID, secondVertexID);
}

void Participant::setMeshEdges(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> vertices)
{
  _impl->setMeshEdges(toSV(meshName), vertices);
}

void Participant::setMeshTriangle(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID)
{
  _impl->setMeshTriangle(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID);
}

void Participant::setMeshTriangles(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> vertices)
{
  _impl->setMeshTriangles(toSV(meshName), vertices);
}

void Participant::setMeshQuad(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID,
    int                    fourthVertexID)
{
  _impl->setMeshQuad(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID,
                     fourthVertexID);
}

void Participant::setMeshQuads(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> vertices)
{
  _impl->setMeshQuads(toSV(meshName), vertices);
}

void Participant::setMeshTetrahedron(
    ::precice::string_view meshName,
    int                    firstVertexID,
    int                    secondVertexID,
    int                    thirdVertexID,
    int                    fourthVertexID)
{
  _impl->setMeshTetrahedron(toSV(meshName), firstVertexID, secondVertexID, thirdVertexID,
                            fourthVertexID);
}

void Participant::setMeshTetrahedra(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> vertices)
{
  _impl->setMeshTetrahedra(toSV(meshName), vertices);
}

void Participant::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   values)
{
  _impl->writeData(toSV(meshName), toSV(dataName), vertices, values);
}

void Participant::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  _impl->readData(toSV(meshName), toSV(dataName), vertices, relativeReadTime, values);
}

void SolverInterface::writeGlobalVectorData(
    std::string_view dataName,
    const double *   value)
{
  _impl->writeGlobalVectorData(dataName, value);
}

void SolverInterface::writeGlobalScalarData(
    std::string_view dataName,
    double           value)
{
  _impl->writeGlobalScalarData(dataName, value);
}

void SolverInterface::readGlobalVectorData(
    std::string_view dataName,
    double *         value) const
{
  _impl->readGlobalVectorData(dataName, value);
}

void SolverInterface::readGlobalVectorData(
    std::string_view dataName,
    double           relativeReadTime,
    double *         value) const
{
  // @todo: needs testing!
  _impl->readGlobalVectorData(dataName, relativeReadTime, value);
}

void SolverInterface::readGlobalScalarData(
    std::string_view dataName,
    double &         value) const
{
  _impl->readGlobalScalarData(dataName, value);
}

void SolverInterface::readGlobalScalarData(
    std::string_view dataName,
    double           relativeReadTime,
    double &         value) const
{
  _impl->readGlobalScalarData(dataName, relativeReadTime, value);
}


void Participant::setMeshAccessRegion(::precice::string_view        meshName,
                                      ::precice::span<const double> boundingBox) const
{
  _impl->setMeshAccessRegion(toSV(meshName), boundingBox);
}

void Participant::getMeshVerticesAndIDs(::precice::string_view    meshName,
                                        ::precice::span<VertexID> ids,
                                        ::precice::span<double>   coordinates) const
{
  _impl->getMeshVerticesAndIDs(toSV(meshName), ids, coordinates);
}

void Participant::writeGradientData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   gradients)
{
  _impl->writeGradientData(toSV(meshName), toSV(dataName), vertices, gradients);
}

} // namespace precice

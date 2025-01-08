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

bool Participant::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(toSV(meshName));
}

void Participant::resetMesh(::precice::string_view meshName)
{
  return _impl->resetMesh(toSV(meshName));
}

bool Participant::requiresGradientDataFor(::precice::string_view meshName,
                                          ::precice::string_view dataName) const
{
  return _impl->requiresGradientDataFor(toSV(meshName), toSV(dataName));
}

VertexID Participant::setMeshVertex(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates)
{
  return _impl->setMeshVertex(toSV(meshName), coordinates);
}

int Participant::getMeshVertexSize(
    ::precice::string_view meshName) const
{
  return _impl->getMeshVertexSize(toSV(meshName));
}

void Participant::setMeshVertices(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates,
    ::precice::span<VertexID>     ids)
{
  _impl->setMeshVertices(toSV(meshName), coordinates, ids);
}

void Participant::setMeshEdge(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second)
{
  _impl->setMeshEdge(toSV(meshName), first, second);
}

void Participant::setMeshEdges(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshEdges(toSV(meshName), ids);
}

void Participant::setMeshTriangle(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third)
{
  _impl->setMeshTriangle(toSV(meshName), first, second, third);
}

void Participant::setMeshTriangles(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshTriangles(toSV(meshName), ids);
}

void Participant::setMeshQuad(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  _impl->setMeshQuad(toSV(meshName), first, second, third, fourth);
}

void Participant::setMeshQuads(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshQuads(toSV(meshName), ids);
}

void Participant::setMeshTetrahedron(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  _impl->setMeshTetrahedron(toSV(meshName), first, second, third, fourth);
}

void Participant::setMeshTetrahedra(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshTetrahedra(toSV(meshName), ids);
}

void Participant::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   values)
{
  _impl->writeData(toSV(meshName), toSV(dataName), ids, values);
}

void Participant::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  _impl->readData(toSV(meshName), toSV(dataName), ids, relativeReadTime, values);
}

void Participant::mapAndReadData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    double                        relativeReadTime,
    ::precice::span<double>       values) const
{
  _impl->mapAndReadData(toSV(meshName), toSV(dataName), coordinates, relativeReadTime, values);
}

void Participant::mapAndWriteData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  _impl->mapAndWriteData(toSV(meshName), toSV(dataName), coordinates, values);
}

void Participant::setMeshAccessRegion(::precice::string_view        meshName,
                                      ::precice::span<const double> boundingBox) const
{
  _impl->setMeshAccessRegion(toSV(meshName), boundingBox);
}

void Participant::getMeshVertexIDsAndCoordinates(::precice::string_view    meshName,
                                                 ::precice::span<VertexID> ids,
                                                 ::precice::span<double>   coordinates) const
{
  _impl->getMeshVertexIDsAndCoordinates(toSV(meshName), ids, coordinates);
}

void Participant::writeGradientData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   gradients)
{
  _impl->writeGradientData(toSV(meshName), toSV(dataName), ids, gradients);
}

} // namespace precice

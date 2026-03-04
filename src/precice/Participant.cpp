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
    void                  *communicator)
    : _impl(new impl::ParticipantImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, {communicator}))
{
}

Participant::~Participant() = default;

void Participant::initialize()
{
  std::scoped_lock lock{_mutex};
  _impl->initialize();
}

void Participant::advance(
    double computedTimeStepSize)
{
  std::scoped_lock lock{_mutex};
  _impl->advance(computedTimeStepSize);
}

void Participant::finalize()
{
  std::scoped_lock lock{_mutex};
  return _impl->finalize();
}

int Participant::getMeshDimensions(::precice::string_view meshName) const
{
  std::scoped_lock lock{_mutex};
  return _impl->getMeshDimensions(toSV(meshName));
}

int Participant::getDataDimensions(::precice::string_view meshName, ::precice::string_view dataName) const
{
  std::scoped_lock lock{_mutex};
  return _impl->getDataDimensions(toSV(meshName), toSV(dataName));
}

bool Participant::isCouplingOngoing() const
{
  std::scoped_lock lock{_mutex};
  return _impl->isCouplingOngoing();
}

bool Participant::isTimeWindowComplete() const
{
  std::scoped_lock lock{_mutex};
  return _impl->isTimeWindowComplete();
}

double Participant::getMaxTimeStepSize() const
{
  std::scoped_lock lock{_mutex};
  return _impl->getMaxTimeStepSize();
}

bool Participant::requiresInitialData()
{
  std::scoped_lock lock{_mutex};
  return _impl->requiresInitialData();
}

bool Participant::requiresReadingCheckpoint()
{
  std::scoped_lock lock{_mutex};
  return _impl->requiresReadingCheckpoint();
}

bool Participant::requiresWritingCheckpoint()
{
  std::scoped_lock lock{_mutex};
  return _impl->requiresWritingCheckpoint();
}

bool Participant::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  std::scoped_lock lock{_mutex};
  return _impl->requiresMeshConnectivityFor(toSV(meshName));
}

void Participant::resetMesh(::precice::string_view meshName)
{
  std::scoped_lock lock{_mutex};
  return _impl->resetMesh(toSV(meshName));
}

bool Participant::requiresGradientDataFor(::precice::string_view meshName,
                                          ::precice::string_view dataName) const
{
  std::scoped_lock lock{_mutex};
  return _impl->requiresGradientDataFor(toSV(meshName), toSV(dataName));
}

VertexID Participant::setMeshVertex(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates)
{
  std::scoped_lock lock{_mutex};
  return _impl->setMeshVertex(toSV(meshName), coordinates);
}

int Participant::getMeshVertexSize(
    ::precice::string_view meshName) const
{
  std::scoped_lock lock{_mutex};
  return _impl->getMeshVertexSize(toSV(meshName));
}

void Participant::setMeshVertices(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates,
    ::precice::span<VertexID>     ids)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshVertices(toSV(meshName), coordinates, ids);
}

void Participant::setMeshEdge(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshEdge(toSV(meshName), first, second);
}

void Participant::setMeshEdges(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshEdges(toSV(meshName), ids);
}

void Participant::setMeshTriangle(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshTriangle(toSV(meshName), first, second, third);
}

void Participant::setMeshTriangles(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshTriangles(toSV(meshName), ids);
}

void Participant::setMeshQuad(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshQuad(toSV(meshName), first, second, third, fourth);
}

void Participant::setMeshQuads(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshQuads(toSV(meshName), ids);
}

void Participant::setMeshTetrahedron(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshTetrahedron(toSV(meshName), first, second, third, fourth);
}

void Participant::setMeshTetrahedra(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshTetrahedra(toSV(meshName), ids);
}

void Participant::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   values)
{
  std::scoped_lock lock{_mutex};
  _impl->writeData(toSV(meshName), toSV(dataName), ids, values);
}

void Participant::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  std::scoped_lock lock{_mutex};
  _impl->readData(toSV(meshName), toSV(dataName), ids, relativeReadTime, values);
}

void Participant::mapAndReadData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    double                        relativeReadTime,
    ::precice::span<double>       values) const
{
  std::scoped_lock lock{_mutex};
  _impl->mapAndReadData(toSV(meshName), toSV(dataName), coordinates, relativeReadTime, values);
}

void Participant::writeAndMapData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  std::scoped_lock lock{_mutex};
  _impl->writeAndMapData(toSV(meshName), toSV(dataName), coordinates, values);
}

void Participant::setMeshAccessRegion(::precice::string_view        meshName,
                                      ::precice::span<const double> boundingBox) const
{
  std::scoped_lock lock{_mutex};
  _impl->setMeshAccessRegion(toSV(meshName), boundingBox);
}

void Participant::getMeshVertexIDsAndCoordinates(::precice::string_view    meshName,
                                                 ::precice::span<VertexID> ids,
                                                 ::precice::span<double>   coordinates) const
{
  std::scoped_lock lock{_mutex};
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

void Participant::startProfilingSection(::precice::string_view sectionName)
{
  _impl->startProfilingSection(toSV(sectionName));
}

void Participant::stopLastProfilingSection()
{
  _impl->stopLastProfilingSection();
}

} // namespace precice

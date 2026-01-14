#include "ParticipantState.hpp"
#include <algorithm>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "WatchIntegral.hpp"
#include "WatchPoint.hpp"
#include "action/Action.hpp"
#include "io/Export.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/String.hpp"
#include "utils/assertion.hpp"
#include "utils/fmt.hpp"

namespace precice::impl {

ParticipantState::ParticipantState(
    std::string                 name,
    mesh::PtrMeshConfiguration &meshConfig)
    : _name(std::move(name))
{
}

/// Configuration interface

void ParticipantState::addAction(action::PtrAction &&action)
{
  auto &context = meshContext(action->getMesh()->getName());
  context.require(action->getMeshRequirement());
  _actions.push_back(std::move(action));
}

void ParticipantState::setUsePrimaryRank(bool useIntraComm)
{
  _useIntraComm = useIntraComm;
}

void ParticipantState::addWatchPoint(
    const PtrWatchPoint &watchPoint)
{
  _watchPoints.push_back(watchPoint);
}

void ParticipantState::addWatchIntegral(
    const PtrWatchIntegral &watchIntegral)
{
  _watchIntegrals.push_back(watchIntegral);
}

void ParticipantState::provideMesh(mesh::PtrMesh mesh)
{
  std::string meshName = mesh->getName();
  PRECICE_TRACE(_name, meshName);
  checkDuplicatedUse(meshName);

  _providedMeshContexts.emplace_back();
  _providedMeshContexts.back().mesh = std::move(mesh);
  _meshContexts[meshName]           = &_providedMeshContexts.back();
  _usedMeshContexts.push_back(&_providedMeshContexts.back());
}

void ParticipantState::receiveMesh(mesh::PtrMesh                                 mesh,
                                   const std::string                            &fromParticipant,
                                   double                                        safetyFactor,
                                   partition::ReceivedPartition::GeometricFilter geoFilter,
                                   const bool                                    allowDirectAccess)
{
  std::string meshName = mesh->getName();
  PRECICE_TRACE(_name, meshName);
  checkDuplicatedUse(meshName);
  PRECICE_ASSERT(!fromParticipant.empty());
  PRECICE_ASSERT(safetyFactor >= 0);

  _receivedMeshContexts.emplace_back();
  auto &context             = _receivedMeshContexts.back();
  context.mesh              = std::move(mesh);
  context.receiveMeshFrom   = fromParticipant;
  context.safetyFactor      = safetyFactor;
  context.geoFilter         = geoFilter;
  context.allowDirectAccess = allowDirectAccess;

  _meshContexts[meshName] = &context;
  _usedMeshContexts.push_back(&context);
}

void ParticipantState::addWriteData(
    const mesh::PtrData &data,
    const mesh::PtrMesh &mesh)
{
  checkDuplicatedData(mesh->getName(), data->getName());
  _writeDataContexts.emplace(MeshDataKey{mesh->getName(), data->getName()}, WriteDataContext(data, mesh));
}

void ParticipantState::addReadData(
    const mesh::PtrData &data,
    const mesh::PtrMesh &mesh)
{
  checkDuplicatedData(mesh->getName(), data->getName());
  _readDataContexts.emplace(MeshDataKey{mesh->getName(), data->getName()}, ReadDataContext(data, mesh));
}

void ParticipantState::addReadMappingContext(
    const MappingContext &mappingContext)
{
  _readMappingContexts.push_back(mappingContext);
}

void ParticipantState::addWriteMappingContext(
    const MappingContext &mappingContext)
{
  _writeMappingContexts.push_back(mappingContext);
}

// Data queries
const ReadDataContext &ParticipantState::readDataContext(std::string_view mesh, std::string_view data) const
{
  auto it = _readDataContexts.find(MeshDataKey{mesh, data});
  PRECICE_CHECK(it != _readDataContexts.end(), "Data \"{}\" does not exist for mesh \"{}\".", data, mesh);
  return it->second;
}

ReadDataContext &ParticipantState::readDataContext(std::string_view mesh, std::string_view data)
{
  auto it = _readDataContexts.find(MeshDataKey{mesh, data});
  PRECICE_CHECK(it != _readDataContexts.end(), "Data \"{}\" does not exist for mesh \"{}\".", data, mesh);
  return it->second;
}

mesh::PtrMesh ParticipantState::findMesh(std::string_view data) const
{
  for (const auto &meshContext : _meshContexts) {
    const auto              &mesh = meshContext.second->mesh->getName();
    MeshDataKey<std::string> key{mesh, std::string{data}};
    const auto               it = _readDataContexts.find(key);
    if (it != _readDataContexts.end()) {
      return meshContext.second->mesh;
    }
  }
  return nullptr;
}

const WriteDataContext &ParticipantState::writeDataContext(std::string_view mesh, std::string_view data) const
{
  auto it = _writeDataContexts.find(MeshDataKey{mesh, data});
  PRECICE_CHECK(it != _writeDataContexts.end(), "Data \"{}\" does not exist in write direction.", data);
  return it->second;
}

WriteDataContext &ParticipantState::writeDataContext(std::string_view mesh, std::string_view data)
{
  auto it = _writeDataContexts.find(MeshDataKey{mesh, data});
  PRECICE_CHECK(it != _writeDataContexts.end(), "Data \"{}\" does not exist in write direction.", data);
  return it->second;
}

bool ParticipantState::hasData(std::string_view mesh, std::string_view data) const
{
  return std::any_of(
      _meshContexts.begin(), _meshContexts.end(),
      [data](const auto &mckv) {
        const auto &meshData = mckv.second->mesh->data();
        return std::any_of(meshData.begin(), meshData.end(), [data](const auto &dptr) {
          return dptr->getName() == data;
        });
      });
}

bool ParticipantState::isDataUsed(std::string_view mesh, std::string_view data) const
{
  const auto &meshData = meshContext(mesh).mesh->data();
  const auto  match    = std::find_if(meshData.begin(), meshData.end(), [data](auto &dptr) { return dptr->getName() == data; });
  return match != meshData.end();
}

bool ParticipantState::isDataRead(std::string_view mesh, std::string_view data) const
{
  return _readDataContexts.count(MeshDataKey{mesh, data}) > 0;
}

bool ParticipantState::isDataWrite(std::string_view mesh, std::string_view data) const
{
  return _writeDataContexts.count(MeshDataKey{mesh, data}) > 0;
}

/// Mesh queries

const MeshContext &ParticipantState::meshContext(std::string_view mesh) const
{
  auto pos = _meshContexts.find(mesh);
  PRECICE_ASSERT(pos != _meshContexts.end());
  return *pos->second;
}

MeshContext &ParticipantState::meshContext(std::string_view mesh)
{
  auto pos = _meshContexts.find(mesh);
  PRECICE_ASSERT(pos != _meshContexts.end());
  return *pos->second;
}

const std::vector<MeshContextVariant> &ParticipantState::usedMeshContexts() const
{
  return _usedMeshContexts;
}

std::vector<MeshContextVariant> &ParticipantState::usedMeshContexts()
{
  return _usedMeshContexts;
}

bool ParticipantState::hasMesh(std::string_view mesh) const
{
  return _meshContexts.count(mesh) > 0;
}

bool ParticipantState::isMeshUsed(std::string_view mesh) const
{
  return std::any_of(
      _usedMeshContexts.begin(), _usedMeshContexts.end(),
      [mesh](const MeshContextVariant &variant) {
        return getMesh(variant).getName() == mesh;
      });
}

bool ParticipantState::isMeshProvided(std::string_view mesh) const
{
  return std::any_of(_providedMeshContexts.begin(), _providedMeshContexts.end(),
                     [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
}

bool ParticipantState::isMeshReceived(std::string_view mesh) const
{
  return std::any_of(_receivedMeshContexts.begin(), _receivedMeshContexts.end(),
                     [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
}

bool ParticipantState::isDirectAccessAllowed(std::string_view mesh) const
{
  auto it = std::find_if(_receivedMeshContexts.begin(), _receivedMeshContexts.end(),
                         [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
  if (it != _receivedMeshContexts.end()) {
    return it->allowDirectAccess;
  }
  return false;
}

ReceivedMeshContext &ParticipantState::receivedMeshContext(std::string_view mesh)
{
  auto it = std::find_if(_receivedMeshContexts.begin(), _receivedMeshContexts.end(),
                         [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
  PRECICE_ASSERT(it != _receivedMeshContexts.end(), "Mesh \"{}\" is not a received mesh", mesh);
  return *it;
}

const ReceivedMeshContext &ParticipantState::receivedMeshContext(std::string_view mesh) const
{
  auto it = std::find_if(_receivedMeshContexts.begin(), _receivedMeshContexts.end(),
                         [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
  PRECICE_ASSERT(it != _receivedMeshContexts.end(), "Mesh \"{}\" is not a received mesh", mesh);
  return *it;
}

ProvidedMeshContext &ParticipantState::providedMeshContext(std::string_view mesh)
{
  auto it = std::find_if(_providedMeshContexts.begin(), _providedMeshContexts.end(),
                         [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
  PRECICE_ASSERT(it != _providedMeshContexts.end(), "Mesh \"{}\" is not a provided mesh", mesh);
  return *it;
}

const ProvidedMeshContext &ParticipantState::providedMeshContext(std::string_view mesh) const
{
  auto it = std::find_if(_providedMeshContexts.begin(), _providedMeshContexts.end(),
                         [mesh](const auto &ctx) { return ctx.mesh->getName() == mesh; });
  PRECICE_ASSERT(it != _providedMeshContexts.end(), "Mesh \"{}\" is not a provided mesh", mesh);
  return *it;
}

std::deque<ProvidedMeshContext> &ParticipantState::providedMeshContexts()
{
  return _providedMeshContexts;
}

const std::deque<ProvidedMeshContext> &ParticipantState::providedMeshContexts() const
{
  return _providedMeshContexts;
}

std::deque<ReceivedMeshContext> &ParticipantState::receivedMeshContexts()
{
  return _receivedMeshContexts;
}

const std::deque<ReceivedMeshContext> &ParticipantState::receivedMeshContexts() const
{
  return _receivedMeshContexts;
}

// Other queries

bool ParticipantState::hasReadMappings() const
{
  return !_readMappingContexts.empty();
}

bool ParticipantState::hasWriteMappings() const
{
  return !_writeMappingContexts.empty();
}

std::vector<MappingContext> &ParticipantState::readMappingContexts()
{
  return _readMappingContexts;
}

std::vector<MappingContext> &ParticipantState::writeMappingContexts()
{
  return _writeMappingContexts;
}

std::vector<action::PtrAction> &ParticipantState::actions()
{
  return _actions;
}

const std::vector<action::PtrAction> &ParticipantState::actions() const
{
  return _actions;
}

void ParticipantState::addExportContext(
    const io::ExportContext &exportContext)
{
  _exportContexts.push_back(exportContext);
}

const std::vector<io::ExportContext> &ParticipantState::exportContexts() const
{
  return _exportContexts;
}

std::vector<PtrWatchPoint> &ParticipantState::watchPoints()
{
  return _watchPoints;
}

std::vector<PtrWatchIntegral> &ParticipantState::watchIntegrals()
{
  return _watchIntegrals;
}

bool ParticipantState::useIntraComm() const
{
  return _useIntraComm;
}

const std::string &ParticipantState::getName() const
{
  return _name;
}

void ParticipantState::exportInitial()
{
  for (const io::ExportContext &context : exportContexts()) {
    if (context.everyNTimeWindows < 1) {
      continue;
    }

    PRECICE_DEBUG("Exporting initial mesh {} to location \"{}\"", context.meshName, context.location);
    context.exporter->doExport(0, 0.0);

    if (context.updateSeries) {
      PRECICE_DEBUG("Exporting series file of mesh {} to location \"{}\"", context.meshName, context.location);
      context.exporter->exportSeries();
    }
  }

  for (const PtrWatchPoint &watchPoint : watchPoints()) {
    watchPoint->exportPointData(0.0);
  }

  for (const PtrWatchIntegral &watchIntegral : watchIntegrals()) {
    watchIntegral->exportIntegralData(0.0);
  }
}

bool ParticipantState::hasExports() const
{
  return !_exportContexts.empty() || !_watchPoints.empty() || !_watchIntegrals.empty();
}

void ParticipantState::exportIntermediate(IntermediateExport exp)
{
  for (const io::ExportContext &context : exportContexts()) {
    if (context.everyIteration) {
      PRECICE_DEBUG("Exporting mesh {} for iteration {} to location \"{}\"", context.meshName, exp.iteration, context.location);
      context.exporter->doExport(exp.iteration, exp.time);
      continue;
    }
    if (exp.complete) {
      PRECICE_DEBUG("Exporting mesh {} for timewindow {} to location \"{}\"", context.meshName, exp.timewindow, context.location);
      context.exporter->doExport(exp.timewindow, exp.time);
    }

    if (exp.final || (exp.complete && context.updateSeries)) {
      PRECICE_DEBUG("Exporting series file of mesh {} to location \"{}\"", context.meshName, context.location);
      context.exporter->exportSeries();
    }
  }

  if (exp.complete) {
    // Export watch point data
    for (const PtrWatchPoint &watchPoint : watchPoints()) {
      watchPoint->exportPointData(exp.time);
    }

    for (const PtrWatchIntegral &watchIntegral : watchIntegrals()) {
      watchIntegral->exportIntegralData(exp.time);
    }
  }
}

// private

void ParticipantState::checkDuplicatedUse(std::string_view mesh)
{
  PRECICE_CHECK(_meshContexts.count(mesh) == 0,
                "Mesh \"{} cannot be used twice by participant {}. "
                "Please remove one of the provide/receive-mesh nodes with name=\"{}\"./>",
                mesh, _name, mesh);
}

void ParticipantState::checkDuplicatedData(std::string_view mesh, std::string_view data)
{
  PRECICE_CHECK(!isDataWrite(mesh, data) && !isDataRead(mesh, data),
                "ParticipantState \"{}\" can read/write data \"{}\" from/to mesh \"{}\" only once. "
                "Please remove any duplicate instances of write-data/read-data nodes.",
                _name, mesh, data);
}

std::string ParticipantState::hintForMesh(std::string_view mesh) const
{
  PRECICE_ASSERT(!hasMesh(mesh));
  PRECICE_ASSERT(!_meshContexts.empty());

  if (_meshContexts.size() == 1) {
    return " This participant only knows mesh \"" + _meshContexts.begin()->first + "\".";
  }

  auto matches = utils::computeMatches(mesh, _meshContexts | boost::adaptors::map_keys);
  if (matches.front().distance < 3) {
    return " Did you mean mesh \"" + matches.front().name + "\"?";
  } else {
    return fmt::format(" Available meshes are: {}", fmt::join(_meshContexts | boost::adaptors::map_keys, ", "));
  }
}

std::string ParticipantState::hintForMeshData(std::string_view mesh, std::string_view data) const
{
  PRECICE_ASSERT(hasMesh(mesh));
  PRECICE_ASSERT(!hasData(mesh, data));
  PRECICE_ASSERT(!_meshContexts.empty());

  // Is there such data in other meshes?
  std::vector<std::string> otherMeshes;
  for (const auto &[_, mc] : _meshContexts) {
    if (mc->mesh->hasDataName(data)) {
      return " Did you mean the data of mesh \"" + mc->mesh->getName() + "\"?";
    }
  }

  // Is there other data in the given mesh?
  auto localData = meshContext(mesh).mesh->availableData();

  if (localData.size() == 1) {
    return " This mesh only knows data \"" + localData.front() + "\".";
  }

  // Was the data typoed?
  auto matches = utils::computeMatches(data, localData);
  if (matches.front().distance < 3) {
    return " Did you mean data \"" + matches.front().name + "\"?";
  }

  return fmt::format(" Available data are: {}", fmt::join(localData, ", "));
}

void ParticipantState::initializeMappingDataCache(std::string_view mappingType)
{
  if (mappingType == "write") {
    for (auto &context : writeDataContexts()) {
      context.initializeMappingDataCache();
    }
  } else {
    for (auto &context : readDataContexts()) {
      context.initializeMappingDataCache();
    }
  }
}

void ParticipantState::configureInputMeshContext(std::string_view fromMesh, impl::MappingContext &mappingContext, mapping::Mapping::MeshRequirement requirement)
{
  meshContext(fromMesh).meshRequirement = std::max(meshContext(fromMesh).meshRequirement, requirement);
  meshContext(fromMesh).fromMappingContexts.push_back(mappingContext);
}

void ParticipantState::configureOutputMeshContext(std::string_view toMesh, impl::MappingContext &mappingContext, mapping::Mapping::MeshRequirement requirement)
{
  meshContext(toMesh).toMappingContexts.push_back(mappingContext);
  meshContext(toMesh).meshRequirement = std::max(meshContext(toMesh).meshRequirement, requirement);
}

} // namespace precice::impl

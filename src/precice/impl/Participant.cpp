#include "Participant.hpp"
#include <algorithm>
#include <ostream>
#include <string>
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
#include "precice/types.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/assertion.hpp"
#include "utils/fmt.hpp"

namespace precice::impl {

Participant::Participant(
    std::string                 name,
    mesh::PtrMeshConfiguration &meshConfig)
    : _name(std::move(name)),
      _meshContexts(meshConfig->meshes().size(), nullptr)
{
}

Participant::~Participant()
{
  for (MeshContext *context : _usedMeshContexts) {
    delete context;
  }
  _usedMeshContexts.clear();
}

/// Configuration interface

void Participant::addAction(action::PtrAction &&action)
{
  auto &context = meshContext(action->getMesh()->getID());
  context.require(action->getMeshRequirement());
  _actions.push_back(std::move(action));
}

void Participant::setUsePrimaryRank(bool useIntraComm)
{
  _useIntraComm = useIntraComm;
}

void Participant::addWatchPoint(
    const PtrWatchPoint &watchPoint)
{
  _watchPoints.push_back(watchPoint);
}

void Participant::addWatchIntegral(
    const PtrWatchIntegral &watchIntegral)
{
  _watchIntegrals.push_back(watchIntegral);
}

void Participant::provideMesh(const mesh::PtrMesh &mesh, bool dynamic)
{
  PRECICE_TRACE(_name, mesh->getName(), mesh->getID());
  checkDuplicatedUse(mesh);

  PRECICE_ASSERT(mesh->getID() < (int) _meshContexts.size());
  auto context         = new MeshContext();
  context->mesh        = mesh;
  context->provideMesh = true;
  if (dynamic) {
    context->dynamic = MeshContext::Dynamicity::Yes;
  }
  _meshContexts[mesh->getID()] = context;
  _usedMeshContexts.push_back(context);
}

void Participant::receiveMesh(const mesh::PtrMesh &                         mesh,
                              const std::string &                           fromParticipant,
                              double                                        safetyFactor,
                              partition::ReceivedPartition::GeometricFilter geoFilter,
                              const bool                                    allowDirectAccess)
{
  PRECICE_TRACE(_name, mesh->getName(), mesh->getID());
  checkDuplicatedUse(mesh);
  PRECICE_ASSERT(mesh->getID() < (int) _meshContexts.size());
  PRECICE_ASSERT(!fromParticipant.empty());
  PRECICE_ASSERT(safetyFactor >= 0);
  auto context               = new MeshContext();
  context->mesh              = mesh;
  context->receiveMeshFrom   = fromParticipant;
  context->safetyFactor      = safetyFactor;
  context->provideMesh       = false;
  context->geoFilter         = geoFilter;
  context->allowDirectAccess = allowDirectAccess;

  _meshContexts[mesh->getID()] = context;

  _usedMeshContexts.push_back(context);
}

void Participant::addWriteData(
    const mesh::PtrData &data,
    const mesh::PtrMesh &mesh)
{
  checkDuplicatedData(data, mesh->getName());
  _writeDataContexts.emplace(data->getID(), WriteDataContext(data, mesh));
}

void Participant::addReadData(
    const mesh::PtrData &data,
    const mesh::PtrMesh &mesh,
    int                  interpolationOrder)
{
  checkDuplicatedData(data, mesh->getName());
  _readDataContexts.emplace(data->getID(), ReadDataContext(data, mesh, interpolationOrder));
}

void Participant::addReadMappingContext(
    const MappingContext &mappingContext)
{
  _readMappingContexts.push_back(mappingContext);
}

void Participant::addWriteMappingContext(
    const MappingContext &mappingContext)
{
  _writeMappingContexts.push_back(mappingContext);
}

// Data queries
const ReadDataContext &Participant::readDataContext(DataID dataID) const
{
  auto it = _readDataContexts.find(dataID);
  PRECICE_CHECK(it != _readDataContexts.end(), "DataID does not exist.")
  return it->second;
}

ReadDataContext &Participant::readDataContext(DataID dataID)
{
  auto it = _readDataContexts.find(dataID);
  PRECICE_CHECK(it != _readDataContexts.end(), "DataID does not exist.")
  return it->second;
}

ReadDataContext &Participant::readDataContext(const std::string &dataName)
{
  auto dataContext = std::find_if(readDataContexts().begin(), readDataContexts().end(), [&dataName](const auto &d) { return d.getDataName() == dataName; });
  PRECICE_ASSERT(dataContext != readDataContexts().end(), "Did not find read data \"{}\".", dataName);

  return *dataContext;
}

const WriteDataContext &Participant::writeDataContext(DataID dataID) const
{
  auto it = _writeDataContexts.find(dataID);
  PRECICE_CHECK(it != _writeDataContexts.end(), "DataID \"{}\" does not exist in write direction.", dataID)
  return it->second;
}

WriteDataContext &Participant::writeDataContext(DataID dataID)
{
  auto it = _writeDataContexts.find(dataID);
  PRECICE_CHECK(it != _writeDataContexts.end(), "DataID \"{}\" does not exist in write direction.", dataID)
  return it->second;
}

bool Participant::hasData(DataID dataID) const
{
  return std::any_of(
      _meshContexts.begin(), _meshContexts.end(),
      [dataID](const auto mcptr) {
        if (!mcptr) {
          return false;
        }
        const auto &meshData = mcptr->mesh->data();
        return std::any_of(meshData.begin(), meshData.end(), [dataID](const auto &dptr) {
          return dptr->getID() == dataID;
        });
      });
}

bool Participant::isDataUsed(const std::string &dataName, MeshID meshID) const
{
  const auto &meshData = meshContext(meshID).mesh->data();
  const auto  match    = std::find_if(meshData.begin(), meshData.end(), [&dataName](auto &dptr) { return dptr->getName() == dataName; });
  return match != meshData.end();
}

bool Participant::isDataRead(DataID dataID) const
{
  return _readDataContexts.count(dataID) > 0;
}

bool Participant::isDataWrite(DataID dataID) const
{
  return _writeDataContexts.count(dataID) > 0;
}

int Participant::getUsedDataID(const std::string &dataName, MeshID meshID) const
{
  const auto &dptr = usedMeshContext(meshID).mesh->data(dataName);
  PRECICE_ASSERT(dptr != nullptr);
  return dptr->getID();
}

std::string Participant::getDataName(DataID dataID) const
{
  for (const MeshContext *mcptr : _meshContexts) {
    if (!mcptr) {
      continue;
    }
    for (const auto &dptr : mcptr->mesh->data()) {
      if (dptr->getID() == dataID) {
        return dptr->getName();
      }
    }
  }
  PRECICE_UNREACHABLE("The dataID {} is invalid.", dataID);
}

/// Mesh queries

const MeshContext &Participant::meshContext(MeshID meshID) const
{
  PRECICE_ASSERT((meshID >= 0) && (meshID < (int) _meshContexts.size()));
  PRECICE_ASSERT(_meshContexts[meshID] != nullptr);
  return *_meshContexts[meshID];
}

MeshContext &Participant::meshContext(MeshID meshID)
{
  PRECICE_TRACE(meshID, _meshContexts.size());
  PRECICE_ASSERT((meshID >= 0) && (meshID < (int) _meshContexts.size()),
                 meshID, _meshContexts.size());
  PRECICE_ASSERT(_meshContexts[meshID] != nullptr);
  return *_meshContexts[meshID];
}

const std::vector<MeshContext *> &Participant::usedMeshContexts() const
{
  return _usedMeshContexts;
}

std::vector<MeshContext *> &Participant::usedMeshContexts()
{
  return _usedMeshContexts;
}

MeshContext &Participant::usedMeshContext(MeshID meshID)
{
  auto pos = std::find_if(_usedMeshContexts.begin(), _usedMeshContexts.end(),
                          [meshID](MeshContext const *context) {
                            return context->mesh->getID() == meshID;
                          });
  PRECICE_ASSERT(pos != _usedMeshContexts.end());
  return **pos;
}

MeshContext const &Participant::usedMeshContext(MeshID meshID) const
{
  auto pos = std::find_if(_usedMeshContexts.begin(), _usedMeshContexts.end(),
                          [meshID](MeshContext const *context) {
                            return context->mesh->getID() == meshID;
                          });
  PRECICE_ASSERT(pos != _usedMeshContexts.end());
  return **pos;
}

MeshContext &Participant::usedMeshContext(const std::string &name)
{
  auto pos = std::find_if(_usedMeshContexts.begin(), _usedMeshContexts.end(),
                          [&name](MeshContext const *context) {
                            return context->mesh->getName() == name;
                          });
  PRECICE_ASSERT(pos != _usedMeshContexts.end());
  return **pos;
}

MeshContext const &Participant::usedMeshContext(const std::string &name) const
{
  auto pos = std::find_if(_usedMeshContexts.begin(), _usedMeshContexts.end(),
                          [&name](MeshContext const *context) {
                            return context->mesh->getName() == name;
                          });
  PRECICE_ASSERT(pos != _usedMeshContexts.end());
  return **pos;
}

bool Participant::hasMesh(MeshID meshID) const
{
  return meshID < static_cast<int>(_meshContexts.size()) && _meshContexts.at(meshID) != nullptr;
}

bool Participant::hasMesh(const std::string &meshName) const
{
  return std::any_of(
      _meshContexts.begin(), _meshContexts.end(),
      [&meshName](const MeshContext *mcptr) {
        return mcptr && meshName == mcptr->mesh->getName();
      });
}

bool Participant::isMeshUsed(MeshID meshID) const
{
  return std::any_of(
      _usedMeshContexts.begin(), _usedMeshContexts.end(),
      [meshID](const MeshContext *mcptr) {
        return mcptr->mesh->getID() == meshID;
      });
}

bool Participant::isMeshUsed(const std::string &meshName) const
{
  return std::any_of(
      _usedMeshContexts.begin(), _usedMeshContexts.end(),
      [&meshName](const MeshContext *mcptr) {
        return mcptr->mesh->getName() == meshName;
      });
}

bool Participant::isMeshProvided(MeshID meshID) const
{
  PRECICE_ASSERT((meshID >= 0) && (meshID < (int) _meshContexts.size()));
  auto context = _meshContexts[meshID];
  return (context != nullptr) && context->provideMesh;
}

bool Participant::isMeshReceived(const std::string &meshName) const
{
  PRECICE_ASSERT(hasMesh(meshName));
  return !usedMeshContext(meshName).provideMesh;
}

bool Participant::isMeshProvided(const std::string &meshName) const
{
  PRECICE_ASSERT(hasMesh(meshName));
  return usedMeshContext(meshName).provideMesh;
}

int Participant::getUsedMeshID(const std::string &meshName) const
{
  return usedMeshContext(meshName).mesh->getID();
}

bool Participant::isDirectAccessAllowed(const int meshID) const
{
  PRECICE_ASSERT((meshID >= 0) && (meshID < (int) _meshContexts.size()));
  auto context = _meshContexts[meshID];
  return context->allowDirectAccess;
}

std::string Participant::getMeshName(MeshID meshID) const
{
  return meshContext(meshID).mesh->getName();
}

std::string Participant::getMeshNameFromData(DataID dataID) const
{
  for (const MeshContext *mcptr : _meshContexts) {
    for (const auto &dptr : mcptr->mesh->data()) {
      if (dptr->getID() == dataID) {
        return mcptr->mesh->getName();
      }
    }
  }
  PRECICE_UNREACHABLE("The dataID {} is invalid.", dataID);
}

// Other queries

std::vector<MappingContext> &Participant::readMappingContexts()
{
  return _readMappingContexts;
}

std::vector<MappingContext> &Participant::writeMappingContexts()
{
  return _writeMappingContexts;
}

std::vector<action::PtrAction> &Participant::actions()
{
  return _actions;
}

const std::vector<action::PtrAction> &Participant::actions() const
{
  return _actions;
}

void Participant::addExportContext(
    const io::ExportContext &exportContext)
{
  _exportContexts.push_back(exportContext);
}

const std::vector<io::ExportContext> &Participant::exportContexts() const
{
  return _exportContexts;
}

bool Participant::isDynamic() const
{
  return std::any_of(
      _usedMeshContexts.begin(), _usedMeshContexts.end(),
      [](auto context) { return context->provideMesh && context->dynamic != MeshContext::Dynamicity::No; });
}

std::set<std::string> Participant::dynamicParticipants() const
{
  std::set<std::string> names;
  if (isDynamic()) {
    names.insert(getName());
  }
  for (const MeshContext *meshContext : usedMeshContexts()) {
    if (meshContext->provideMesh || meshContext->dynamic == MeshContext::Dynamicity::No) {
      continue;
    }
    names.insert(meshContext->receiveMeshFrom);
  }
  PRECICE_WARN("Dynamic Participants are {}", names);
  return names;
}

std::vector<PtrWatchPoint> &Participant::watchPoints()
{
  return _watchPoints;
}

std::vector<PtrWatchIntegral> &Participant::watchIntegrals()
{
  return _watchIntegrals;
}

bool Participant::useIntraComm() const
{
  return _useIntraComm;
}

const std::string &Participant::getName() const
{
  return _name;
}

void Participant::exportInitial()
{
  for (const io::ExportContext &context : exportContexts()) {
    if (context.everyNTimeWindows < 1) {
      continue;
    }

    for (const MeshContext *meshContext : usedMeshContexts()) {
      auto &mesh = *meshContext->mesh;
      PRECICE_DEBUG("Exporting initial mesh {} to location \"{}\"", mesh.getName(), context.location);
      context.exporter->doExport(fmt::format("{}-{}.init", mesh.getName(), getName()), context.location, mesh);
    }
  }
}

void Participant::exportFinal()
{
  for (const io::ExportContext &context : exportContexts()) {
    if (context.everyNTimeWindows < 1) {
      continue;
    }

    for (const MeshContext *meshContext : usedMeshContexts()) {
      auto &mesh = *meshContext->mesh;
      PRECICE_DEBUG("Exporting final mesh {} to location \"{}\"", mesh.getName(), context.location);
      context.exporter->doExport(fmt::format("{}-{}.final", mesh.getName(), getName()), context.location, mesh);
    }
  }
}

void Participant::exportIntermediate(IntermediateExport exp)
{
  for (const io::ExportContext &context : exportContexts()) {
    if (exp.complete && (context.everyNTimeWindows > 0) && (exp.timewindow % context.everyNTimeWindows == 0)) {
      for (const MeshContext *meshContext : usedMeshContexts()) {
        auto &mesh = *meshContext->mesh;
        PRECICE_DEBUG("Exporting mesh {} for timewindow {} to location \"{}\"", mesh.getName(), exp.timewindow, context.location);
        context.exporter->doExport(fmt::format("{}-{}.dt{}", mesh.getName(), getName(), exp.timewindow), context.location, mesh);
      }
    }

    if (context.everyIteration) {
      for (const MeshContext *meshContext : usedMeshContexts()) {
        auto &mesh = *meshContext->mesh;
        PRECICE_DEBUG("Exporting mesh {} for iteration {} to location \"{}\"", meshContext->mesh->getName(), exp.iteration, context.location);
        /// @todo this is the global iteration count. Shouldn't this be local to the timestep? example .dtN.itM or similar
        context.exporter->doExport(fmt::format("{}-{}.it{}", mesh.getName(), getName(), exp.iteration), context.location, mesh);
      }
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

void Participant::checkDuplicatedUse(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT((int) _meshContexts.size() > mesh->getID());
  PRECICE_CHECK(_meshContexts[mesh->getID()] == nullptr,
                "Mesh \"{} cannot be used twice by participant {}. "
                "Please remove one of the provide/receive-mesh nodes with name=\"{}\"./>",
                mesh->getName(), _name, mesh->getName());
}

void Participant::checkDuplicatedData(const mesh::PtrData &data, const std::string &meshName)
{
  PRECICE_CHECK(!isDataWrite(data->getID()) && !isDataRead(data->getID()),
                "Participant \"{}\" can read/write data \"{}\" from/to mesh \"{}\" only once. "
                "Please remove any duplicate instances of write-data/read-data nodes.",
                _name, meshName, data->getName());
}

} // namespace precice::impl

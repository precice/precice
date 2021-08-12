#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <functional>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream>
#include <tuple>
#include <utility>

#include "SolverInterfaceImpl.hpp"
#include "action/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "io/Export.hpp"
#include "io/ExportContext.hpp"
#include "io/SharedPointer.hpp"
#include "logging/LogConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "partition/Partition.hpp"
#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"
#include "partition/SharedPointer.hpp"
#include "precice/config/Configuration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "precice/config/SharedPointer.hpp"
#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/impl/CommonErrorMessages.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/MappingContext.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/ValidationMacros.hpp"
#include "precice/impl/WatchIntegral.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "precice/impl/versions.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/EigenIO.hpp"
#include "utils/Event.hpp"
#include "utils/EventUtils.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/PointerVector.hpp"
#include "utils/algorithm.hpp"
#include "utils/assertion.hpp"
#include "xml/XMLTag.hpp"

using precice::utils::Event;
using precice::utils::EventRegistry;

namespace precice {

/// Enabled further inter- and intra-solver synchronisation
bool syncMode = false;

namespace impl {

SolverInterfaceImpl::SolverInterfaceImpl(
    std::string        participantName,
    const std::string &configurationFileName,
    int                accessorProcessRank,
    int                accessorCommunicatorSize,
    void *             communicator)
    : _accessorName(std::move(participantName)),
      _accessorProcessRank(accessorProcessRank),
      _accessorCommunicatorSize(accessorCommunicatorSize)
{
  PRECICE_CHECK(!_accessorName.empty(),
                "This participant's name is an empty string. "
                "When constructing a preCICE interface you need to pass the name of the "
                "participant as first argument to the constructor.");
  PRECICE_CHECK(_accessorProcessRank >= 0,
                "The solver process index needs to be a non-negative number, not: {}. "
                "Please check the value given when constructing a preCICE interface.",
                _accessorProcessRank);
  PRECICE_CHECK(_accessorCommunicatorSize >= 1,
                "The solver process size needs to be a positive number, not: {}. "
                "Please check the value given when constructing a preCICE interface.",
                _accessorCommunicatorSize);
  PRECICE_CHECK(_accessorProcessRank < _accessorCommunicatorSize,
                "The solver process index, currently: {}  needs to be smaller than the solver process size, currently: {}. "
                "Please check the values given when constructing a preCICE interface.",
                _accessorProcessRank, _accessorCommunicatorSize);

// Set the global communicator to the passed communicator.
// This is a noop if preCICE is not configured with MPI.
// nullpointer signals to use MPI_COMM_WORLD
#ifndef PRECICE_NO_MPI
  if (communicator != nullptr) {
    auto commptr = static_cast<utils::Parallel::Communicator *>(communicator);
    utils::Parallel::registerUserProvidedComm(*commptr);
  }
#endif

  logging::setParticipant(_accessorName);

  configure(configurationFileName);

// This block cannot be merge with the one above as only configure calls
// utils::Parallel::initializeMPI, which is needed for getProcessRank.
#ifndef PRECICE_NO_MPI
  if (communicator != nullptr) {
    const auto currentRank = utils::Parallel::current()->rank();
    PRECICE_CHECK(_accessorProcessRank == currentRank,
                  "The solver process index given in the preCICE interface constructor({}) does not match the rank of the passed MPI communicator ({}).",
                  _accessorProcessRank, currentRank);
    const auto currentSize = utils::Parallel::current()->size();
    PRECICE_CHECK(_accessorCommunicatorSize == currentSize,
                  "The solver process size given in the preCICE interface constructor({}) does not match the size of the passed MPI communicator ({}).",
                  _accessorCommunicatorSize, currentSize);
  }
#endif
}

SolverInterfaceImpl::SolverInterfaceImpl(
    std::string        participantName,
    const std::string &configurationFileName,
    int                accessorProcessRank,
    int                accessorCommunicatorSize)
    : SolverInterfaceImpl::SolverInterfaceImpl(std::move(participantName), configurationFileName, accessorProcessRank, accessorCommunicatorSize, nullptr)
{
}

SolverInterfaceImpl::~SolverInterfaceImpl()
{
  if (_state != State::Finalized) {
    PRECICE_INFO("Implicitly finalizing in destructor");
    finalize();
  }
}

void SolverInterfaceImpl::configure(
    const std::string &configurationFileName)
{
  config::Configuration config;
  utils::Parallel::initializeManagedMPI(nullptr, nullptr);
  logging::setMPIRank(utils::Parallel::current()->rank());
  xml::ConfigurationContext context{
      _accessorName,
      _accessorProcessRank,
      _accessorCommunicatorSize};
  xml::configure(config.getXMLTag(), context, configurationFileName);
  if (_accessorProcessRank == 0) {
    PRECICE_INFO("This is preCICE version {}", PRECICE_VERSION);
    PRECICE_INFO("Revision info: {}", precice::preciceRevision);
#ifndef NDEBUG
    PRECICE_INFO("Configuration: Debug");
#else
    PRECICE_INFO("Configuration: Release (Debug and Trace log unavailable)");
#endif
    PRECICE_INFO("Configuring preCICE with configuration \"{}\"", configurationFileName);
    PRECICE_INFO("I am participant \"{}\"", _accessorName);
  }
  configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceImpl::configure(
    const config::SolverInterfaceConfiguration &config)
{
  PRECICE_TRACE();

  Event                    e("configure"); // no precice::syncMode as this is not yet configured here
  utils::ScopedEventPrefix sep("configure/");

  mesh::Data::resetDataCount();
  _meshLock.clear();

  _dimensions = config.getDimensions();
  _accessor   = determineAccessingParticipant(config);
  _accessor->setMeshIdManager(config.getMeshConfiguration()->extractMeshIdManager());

  PRECICE_ASSERT(_accessorCommunicatorSize == 1 || _accessor->useMaster(),
                 "A parallel participant needs a master communication");
  PRECICE_CHECK(not(_accessorCommunicatorSize == 1 && _accessor->useMaster()),
                "You cannot use a master communication with a serial participant. "
                "If you do not know exactly what a master communication is and why you want to use it "
                "you probably just want to remove the master tag from the preCICE configuration.");

  utils::MasterSlave::configure(_accessorProcessRank, _accessorCommunicatorSize);

  _participants = config.getParticipantConfiguration()->getParticipants();
  configureM2Ns(config.getM2NConfiguration());

  PRECICE_CHECK(_participants.size() > 1,
                "In the preCICE configuration, only one participant is defined. "
                "One participant makes no coupled simulation. "
                "Please add at least another one.");
  configurePartitions(config.getM2NConfiguration());

  cplscheme::PtrCouplingSchemeConfiguration cplSchemeConfig =
      config.getCouplingSchemeConfiguration();
  _couplingScheme = cplSchemeConfig->getCouplingScheme(_accessorName);

  // Register all MeshIds to the lock, but unlock them straight away as
  // writing is allowed after configuration.
  for (const MeshContext *meshContext : _accessor->usedMeshContexts()) {
    _meshLock.add(meshContext->mesh->getID(), false);
  }

  utils::EventRegistry::instance().initialize("precice-" + _accessorName, "", utils::Parallel::current()->comm);

  PRECICE_DEBUG("Initialize master-slave communication");
  if (utils::MasterSlave::isParallel()) {
    initializeMasterSlaveCommunication();
  }

  auto &solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.start(precice::syncMode);
}

double SolverInterfaceImpl::initialize()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Finalized, "initialize() cannot be called after finalize().")
  PRECICE_CHECK(_state != State::Initialized, "initialize() may only be called once.");
  PRECICE_ASSERT(not _couplingScheme->isInitialized());
  auto &solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.pause(precice::syncMode);
  Event                    e("initialize", precice::syncMode);
  utils::ScopedEventPrefix sep("initialize/");

  // Setup communication

  PRECICE_INFO("Setting up master communication to coupling partner/s");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n       = m2nPair.second;
    bool  requesting = bm2n.isRequesting;
    if (bm2n.m2n->isConnected()) {
      PRECICE_DEBUG("Master connection {} {} already connected.", (requesting ? "from" : "to"), bm2n.remoteName);
    } else {
      PRECICE_DEBUG((requesting ? "Awaiting master connection from {}" : "Establishing master connection to {}"), bm2n.remoteName);
      bm2n.prepareEstablishment();
      bm2n.connectMasters();
      PRECICE_DEBUG("Established master connection {} {}", (requesting ? "from " : "to "), bm2n.remoteName);
    }
  }

  PRECICE_INFO("Masters are connected");

  compareBoundingBoxes();

  PRECICE_INFO("Setting up preliminary slaves communication to coupling partner/s");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n = m2nPair.second;
    bm2n.preConnectSlaves();
  }

  computePartitions();

  PRECICE_INFO("Setting up slaves communication to coupling partner/s");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n = m2nPair.second;
    bm2n.connectSlaves();
    PRECICE_DEBUG("Established slaves connection {} {}", (bm2n.isRequesting ? "from " : "to "), bm2n.remoteName);
  }
  PRECICE_INFO("Slaves are connected");

  for (auto &m2nPair : _m2ns) {
    m2nPair.second.cleanupEstablishment();
  }

  PRECICE_DEBUG("Initialize watchpoints");
  for (PtrWatchPoint &watchPoint : _accessor->watchPoints()) {
    watchPoint->initialize();
  }
  for (PtrWatchIntegral &watchIntegral : _accessor->watchIntegrals()) {
    watchIntegral->initialize();
  }

  // Initialize coupling state, overwrite these values for restart
  double time       = 0.0;
  int    timeWindow = 1;

  PRECICE_DEBUG("Initialize coupling schemes");
  _couplingScheme->initialize(time, timeWindow);
  PRECICE_ASSERT(_couplingScheme->isInitialized());

  double dt = 0.0;

  dt = _couplingScheme->getNextTimestepMaxLength();

  if (_couplingScheme->hasDataBeenReceived()) {
    performDataActions({action::Action::READ_MAPPING_PRIOR}, 0.0, 0.0, 0.0, dt);
    mapReadData();
    performDataActions({action::Action::READ_MAPPING_POST}, 0.0, 0.0, 0.0, dt);
  }

  PRECICE_INFO(_couplingScheme->printCouplingState());

  solverInitEvent.start(precice::syncMode);

  _meshLock.lockAll();

  _state = State::Initialized;

  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl::initializeData()
{
  PRECICE_TRACE();
  PRECICE_CHECK(!_hasInitializedData, "initializeData() may only be called once.");
  PRECICE_CHECK(_state != State::Finalized, "initializeData() cannot be called after finalize().")
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before initializeData()");
  PRECICE_ASSERT(_couplingScheme->isInitialized());
  PRECICE_CHECK(not(_couplingScheme->sendsInitializedData() && isActionRequired(constants::actionWriteInitialData())),
                "Initial data has to be written to preCICE by calling an appropriate write...Data() function before calling initializeData(). "
                "Did you forget to call markActionFulfilled(precice::constants::actionWriteInitialData()) after writing initial data?");

  auto &solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.pause(precice::syncMode);

  Event                    e("initializeData", precice::syncMode);
  utils::ScopedEventPrefix sep("initializeData/");

  PRECICE_DEBUG("Initialize data");
  double dt = _couplingScheme->getNextTimestepMaxLength();

  performDataActions({action::Action::WRITE_MAPPING_PRIOR}, 0.0, 0.0, 0.0, dt);
  mapWrittenData();
  performDataActions({action::Action::WRITE_MAPPING_POST}, 0.0, 0.0, 0.0, dt);

  _couplingScheme->initializeData();

  if (_couplingScheme->hasDataBeenReceived()) {
    performDataActions({action::Action::READ_MAPPING_PRIOR}, 0.0, 0.0, 0.0, dt);
    mapReadData();
    performDataActions({action::Action::READ_MAPPING_POST}, 0.0, 0.0, 0.0, dt);
  }
  resetWrittenData();
  PRECICE_DEBUG("Plot output");
  for (const io::ExportContext &context : _accessor->exportContexts()) {
    if (context.everyNTimeWindows != -1) {
      std::ostringstream suffix;
      suffix << _accessorName << ".init";
      exportMesh(suffix.str());
    }
  }
  solverInitEvent.start(precice::syncMode);

  _hasInitializedData = true;
}

double SolverInterfaceImpl::advance(
    double computedTimestepLength)
{

  PRECICE_TRACE(computedTimestepLength);

  // Events for the solver time, stopped when we enter, restarted when we leave advance
  auto &solverEvent = EventRegistry::instance().getStoredEvent("solver.advance");
  solverEvent.stop(precice::syncMode);
  auto &solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.stop(precice::syncMode);

  Event                    e("advance", precice::syncMode);
  utils::ScopedEventPrefix sep("advance/");

  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before advance().");
  PRECICE_CHECK(_state != State::Finalized, "advance() cannot be called after finalize().")
  PRECICE_ASSERT(_couplingScheme->isInitialized());
  PRECICE_CHECK(isCouplingOngoing(), "advance() cannot be called when isCouplingOngoing() returns false.");
  PRECICE_CHECK((not _couplingScheme->receivesInitializedData() && not _couplingScheme->sendsInitializedData()) || (_hasInitializedData),
                "initializeData() needs to be called before advance if data has to be initialized.");
  PRECICE_CHECK(!math::equals(computedTimestepLength, 0.0), "advance() cannot be called with a timestep size of 0.");
  PRECICE_CHECK(computedTimestepLength > 0.0, "advance() cannot be called with a negative timestep size {}.", computedTimestepLength);
  _numberAdvanceCalls++;

#ifndef NDEBUG
  PRECICE_DEBUG("Synchronize timestep length");
  if (utils::MasterSlave::isParallel()) {
    syncTimestep(computedTimestepLength);
  }
#endif

  double timeWindowSize         = 0.0; // Length of (full) current time window
  double timeWindowComputedPart = 0.0; // Length of computed part of (full) current time window
  double time                   = 0.0; // Current time

  // Update the coupling scheme time state. Necessary to get correct remainder.
  _couplingScheme->addComputedTime(computedTimestepLength);

  if (_couplingScheme->hasTimeWindowSize()) {
    timeWindowSize = _couplingScheme->getTimeWindowSize();
  } else {
    timeWindowSize = computedTimestepLength;
  }
  timeWindowComputedPart = timeWindowSize - _couplingScheme->getThisTimeWindowRemainder();
  time                   = _couplingScheme->getTime();

  if (_couplingScheme->willDataBeExchanged(0.0)) {
    performDataActions({action::Action::WRITE_MAPPING_PRIOR}, time, computedTimestepLength, timeWindowComputedPart, timeWindowSize);
    mapWrittenData();
    performDataActions({action::Action::WRITE_MAPPING_POST}, time, computedTimestepLength, timeWindowComputedPart, timeWindowSize);
  }

  PRECICE_DEBUG("Advance coupling scheme");
  _couplingScheme->advance();

  if (_couplingScheme->hasDataBeenReceived()) {
    performDataActions({action::Action::READ_MAPPING_PRIOR}, time, computedTimestepLength, timeWindowComputedPart, timeWindowSize);
    mapReadData();
    performDataActions({action::Action::READ_MAPPING_POST}, time, computedTimestepLength, timeWindowComputedPart, timeWindowSize);
  }

  if (_couplingScheme->isTimeWindowComplete()) {
    performDataActions({action::Action::ON_TIME_WINDOW_COMPLETE_POST}, time, computedTimestepLength, timeWindowComputedPart, timeWindowSize);
  }

  PRECICE_INFO(_couplingScheme->printCouplingState());

  PRECICE_DEBUG("Handle exports");
  handleExports();

  resetWrittenData();

  _meshLock.lockAll();
  solverEvent.start(precice::syncMode);
  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl::finalize()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Finalized, "finalize() may only be called once.")

  // Events for the solver time, finally stopped here
  auto &solverEvent = EventRegistry::instance().getStoredEvent("solver.advance");
  solverEvent.stop(precice::syncMode);

  Event                    e("finalize"); // no precice::syncMode here as MPI is already finalized at destruction of this event
  utils::ScopedEventPrefix sep("finalize/");

  if (_state == State::Initialized) {

    PRECICE_ASSERT(_couplingScheme->isInitialized());
    PRECICE_DEBUG("Finalize coupling scheme");
    _couplingScheme->finalize();

    PRECICE_DEBUG("Handle exports");
    for (const io::ExportContext &context : _accessor->exportContexts()) {
      if (context.everyNTimeWindows != -1) {
        std::ostringstream suffix;
        suffix << _accessorName << ".final";
        exportMesh(suffix.str());
      }
    }
    closeCommunicationChannels(CloseChannels::All);
  }

  // Release ownership
  _couplingScheme.reset();
  _participants.clear();
  _accessor.reset();

  // Close Connections
  PRECICE_DEBUG("Close master-slave communication");
  if (utils::MasterSlave::isParallel()) {
    utils::MasterSlave::_communication->closeConnection();
    utils::MasterSlave::_communication = nullptr;
  }
  _m2ns.clear();

  // Stop and print Event logging
  e.stop();

  // Finalize PETSc and Events first
  utils::Petsc::finalize();
  utils::EventRegistry::instance().finalize();

  // Printing requires finalization
  if (not precice::utils::MasterSlave::isSlave()) {
    utils::EventRegistry::instance().printAll();
  }

  // Finally clear events and finalize MPI
  utils::EventRegistry::instance().clear();
  utils::Parallel::finalizeManagedMPI();
  _state = State::Finalized;
}

int SolverInterfaceImpl::getDimensions() const
{
  PRECICE_TRACE(_dimensions);
  return _dimensions;
}

bool SolverInterfaceImpl::isCouplingOngoing() const
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isCouplingOngoing() can be evaluated.");
  PRECICE_CHECK(_state != State::Finalized, "isCouplingOngoing() cannot be called after finalize().");
  return _couplingScheme->isCouplingOngoing();
}

bool SolverInterfaceImpl::isReadDataAvailable() const
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isReadDataAvailable().");
  PRECICE_CHECK(_state != State::Finalized, "isReadDataAvailable() cannot be called after finalize().");
  return _couplingScheme->hasDataBeenReceived();
}

bool SolverInterfaceImpl::isWriteDataRequired(
    double computedTimestepLength) const
{
  PRECICE_TRACE(computedTimestepLength);
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isWriteDataRequired().");
  PRECICE_CHECK(_state != State::Finalized, "isWriteDataRequired() cannot be called after finalize().");
  return _couplingScheme->willDataBeExchanged(computedTimestepLength);
}

bool SolverInterfaceImpl::isTimeWindowComplete() const
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isTimeWindowComplete().");
  PRECICE_CHECK(_state != State::Finalized, "isTimeWindowComplete() cannot be called after finalize().");
  return _couplingScheme->isTimeWindowComplete();
}

bool SolverInterfaceImpl::isActionRequired(
    const std::string &action) const
{
  PRECICE_TRACE(action, _couplingScheme->isActionRequired(action));
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isActionRequired(...).");
  PRECICE_CHECK(_state != State::Finalized, "isActionRequired(...) cannot be called after finalize().");
  return _couplingScheme->isActionRequired(action);
}

void SolverInterfaceImpl::markActionFulfilled(
    const std::string &action)
{
  PRECICE_TRACE(action);
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before markActionFulfilled(...).");
  PRECICE_CHECK(_state != State::Finalized, "markActionFulfilled(...) cannot be called after finalize().");
  _couplingScheme->markActionFulfilled(action);
}

bool SolverInterfaceImpl::hasToEvaluateSurrogateModel() const
{
  return false;
}

bool SolverInterfaceImpl::hasToEvaluateFineModel() const
{
  return true;
}

bool SolverInterfaceImpl::hasMesh(
    const std::string &meshName) const
{
  PRECICE_TRACE(meshName);
  return _accessor->hasMesh(meshName);
}

int SolverInterfaceImpl::getMeshID(
    const std::string &meshName) const
{
  PRECICE_TRACE(meshName);
  PRECICE_CHECK(_accessor->hasMesh(meshName),
                "The given mesh name \"{}\" is unknown to preCICE. "
                "Please check the mesh definitions in the configuration.",
                meshName);
  PRECICE_CHECK(_accessor->isMeshUsed(meshName),
                "The given mesh name \"{0}\" is not used by the participant \"{1}\". "
                "Please define a <use-mesh name=\"{0}\"/> node for the particpant \"{1}\".",
                meshName, _accessorName);
  return _accessor->getUsedMeshID(meshName);
}

std::set<int> SolverInterfaceImpl::getMeshIDs() const
{
  PRECICE_TRACE();
  std::set<int> ids;
  for (const impl::MeshContext *context : _accessor->usedMeshContexts()) {
    ids.insert(context->mesh->getID());
  }
  return ids;
}

bool SolverInterfaceImpl::hasData(
    const std::string &dataName, MeshID meshID) const
{
  PRECICE_TRACE(dataName, meshID);
  PRECICE_VALIDATE_MESH_ID(meshID);
  return _accessor->isDataUsed(dataName, meshID);
}

int SolverInterfaceImpl::getDataID(
    const std::string &dataName, MeshID meshID) const
{
  PRECICE_TRACE(dataName, meshID);
  PRECICE_VALIDATE_MESH_ID(meshID);
  PRECICE_CHECK(_accessor->isDataUsed(dataName, meshID),
                "Data with name \"{0}\" is not defined on mesh \"{1}\". "
                "Please add <use-data name=\"{0}\"/> under <mesh name=\"{1}\"/>.",
                dataName, _accessor->getMeshName(meshID));
  return _accessor->getUsedDataID(dataName, meshID);
}

int SolverInterfaceImpl::getMeshVertexSize(
    MeshID meshID) const
{
  PRECICE_TRACE(meshID);
  PRECICE_REQUIRE_MESH_USE(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  PRECICE_ASSERT(context.mesh.get() != nullptr);
  return context.mesh->vertices().size();
}

/// @todo Currently not supported as we would need to re-compute the re-partition
void SolverInterfaceImpl::resetMesh(
    MeshID meshID)
{
  PRECICE_TRACE(meshID);
  PRECICE_VALIDATE_MESH_ID(meshID);
  impl::MeshContext &context = _accessor->usedMeshContext(meshID);
  /*
  bool               hasMapping = context.fromMappingContext.mapping || context.toMappingContext.mapping;
  bool               isStationary =
      context.fromMappingContext.timing == mapping::MappingConfiguration::INITIAL &&
      context.toMappingContext.timing == mapping::MappingConfiguration::INITIAL;
  */

  PRECICE_DEBUG("Clear mesh positions for mesh \"{}\"", context.mesh->getName());
  _meshLock.unlock(meshID);
  context.mesh->clear();
}

int SolverInterfaceImpl::setMeshVertex(
    int           meshID,
    const double *position)
{
  PRECICE_TRACE(meshID);
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  Eigen::VectorXd internalPosition{
      Eigen::Map<const Eigen::VectorXd>{position, _dimensions}};
  PRECICE_DEBUG("Position = {}", internalPosition.format(utils::eigenio::debug()));
  int           index   = -1;
  MeshContext & context = _accessor->usedMeshContext(meshID);
  mesh::PtrMesh mesh(context.mesh);
  PRECICE_DEBUG("MeshRequirement: {}", context.meshRequirement);
  index = mesh->createVertex(internalPosition).getID();
  mesh->allocateDataValues();
  return index;
}

void SolverInterfaceImpl::setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_TRACE(meshID, size);
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext & context = _accessor->usedMeshContext(meshID);
  mesh::PtrMesh mesh(context.mesh);
  PRECICE_DEBUG("Set positions");
  const Eigen::Map<const Eigen::MatrixXd> posMatrix{
      positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
  for (int i = 0; i < size; ++i) {
    Eigen::VectorXd current(posMatrix.col(i));
    ids[i] = mesh->createVertex(current).getID();
  }
  mesh->allocateDataValues();
}

void SolverInterfaceImpl::getMeshVertices(
    int        meshID,
    size_t     size,
    const int *ids,
    double *   positions) const
{
  PRECICE_TRACE(meshID, size);
  PRECICE_REQUIRE_MESH_USE(meshID);
  MeshContext & context = _accessor->usedMeshContext(meshID);
  mesh::PtrMesh mesh(context.mesh);
  PRECICE_DEBUG("Get positions");
  auto &vertices = mesh->vertices();
  PRECICE_ASSERT(size <= vertices.size(), size, vertices.size());
  Eigen::Map<Eigen::MatrixXd> posMatrix{
      positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
  for (size_t i = 0; i < size; i++) {
    const size_t id = ids[i];
    PRECICE_ASSERT(id < vertices.size(), id, vertices.size());
    posMatrix.col(i) = vertices[id].getCoords();
  }
}

void SolverInterfaceImpl::getMeshVertexIDsFromPositions(
    int           meshID,
    size_t        size,
    const double *positions,
    int *         ids) const
{
  PRECICE_TRACE(meshID, size);
  PRECICE_REQUIRE_MESH_USE(meshID);
  MeshContext & context = _accessor->usedMeshContext(meshID);
  mesh::PtrMesh mesh(context.mesh);
  PRECICE_DEBUG("Get IDs");
  const auto &                      vertices = mesh->vertices();
  Eigen::Map<const Eigen::MatrixXd> posMatrix{
      positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
  const auto vsize = vertices.size();
  for (size_t i = 0; i < size; i++) {
    size_t j = 0;
    for (; j < vsize; j++) {
      if (math::equals(posMatrix.col(i), vertices[j].getCoords())) {
        break;
      }
    }
    if (j == vsize) {
      std::ostringstream err;
      err << "Unable to find a vertex on mesh \"" << mesh->getName() << "\" at position (";
      err << posMatrix.col(i)[0] << ", " << posMatrix.col(i)[1];
      if (_dimensions == 3) {
        err << ", " << posMatrix.col(i)[2];
      }
      err << "). The request failed for query " << i + 1 << " out of " << size << '.';
      PRECICE_ERROR(err.str());
    }
    ids[i] = j;
  }
}

int SolverInterfaceImpl::setMeshEdge(
    MeshID meshID,
    int    firstVertexID,
    int    secondVertexID)
{
  PRECICE_TRACE(meshID, firstVertexID, secondVertexID);
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh->isValidVertexID(firstVertexID), errorInvalidVertexID(firstVertexID));
    PRECICE_CHECK(mesh->isValidVertexID(secondVertexID), errorInvalidVertexID(secondVertexID));
    mesh::Vertex &v0 = mesh->vertices()[firstVertexID];
    mesh::Vertex &v1 = mesh->vertices()[secondVertexID];
    return mesh->createEdge(v0, v1).getID();
  }
  return -1;
}

void SolverInterfaceImpl::setMeshTriangle(
    MeshID meshID,
    int    firstEdgeID,
    int    secondEdgeID,
    int    thirdEdgeID)
{
  PRECICE_TRACE(meshID, firstEdgeID,
                secondEdgeID, thirdEdgeID);
  PRECICE_CHECK(_dimensions == 3, "setMeshTriangle is only possible for 3D cases."
                                  " Please set the dimension to 3 in the preCICE configuration file.");
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidEdgeID;
    PRECICE_CHECK(mesh->isValidEdgeID(firstEdgeID), errorInvalidEdgeID(firstEdgeID));
    PRECICE_CHECK(mesh->isValidEdgeID(secondEdgeID), errorInvalidEdgeID(secondEdgeID));
    PRECICE_CHECK(mesh->isValidEdgeID(thirdEdgeID), errorInvalidEdgeID(thirdEdgeID));
    PRECICE_CHECK(utils::unique_elements(utils::make_array(firstEdgeID, secondEdgeID, thirdEdgeID)),
                  "setMeshTriangle() was called with repeated Edge IDs ({}, {}, {}).",
                  firstEdgeID, secondEdgeID, thirdEdgeID);
    mesh::Edge &e0 = mesh->edges()[firstEdgeID];
    mesh::Edge &e1 = mesh->edges()[secondEdgeID];
    mesh::Edge &e2 = mesh->edges()[thirdEdgeID];
    PRECICE_CHECK(e0.connectedTo(e1) && e1.connectedTo(e2) && e2.connectedTo(e0),
                  "setMeshTriangle() was called with Edge IDs ({}, {}, {}), which identify unconnected Edges.",
                  firstEdgeID, secondEdgeID, thirdEdgeID);
    mesh->createTriangle(e0, e1, e2);
  }
}

void SolverInterfaceImpl::setMeshTriangleWithEdges(
    MeshID meshID,
    int    firstVertexID,
    int    secondVertexID,
    int    thirdVertexID)
{
  PRECICE_TRACE(meshID, firstVertexID,
                secondVertexID, thirdVertexID);
  PRECICE_CHECK(_dimensions == 3, "setMeshTriangleWithEdges is only possible for 3D cases."
                                  " Please set the dimension to 3 in the preCICE configuration file.");
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh->isValidVertexID(firstVertexID), errorInvalidVertexID(firstVertexID));
    PRECICE_CHECK(mesh->isValidVertexID(secondVertexID), errorInvalidVertexID(secondVertexID));
    PRECICE_CHECK(mesh->isValidVertexID(thirdVertexID), errorInvalidVertexID(thirdVertexID));
    PRECICE_CHECK(utils::unique_elements(utils::make_array(firstVertexID, secondVertexID, thirdVertexID)),
                  "setMeshTriangleWithEdges() was called with repeated Vertex IDs ({}, {}, {}).",
                  firstVertexID, secondVertexID, thirdVertexID);
    mesh::Vertex *vertices[3];
    vertices[0] = &mesh->vertices()[firstVertexID];
    vertices[1] = &mesh->vertices()[secondVertexID];
    vertices[2] = &mesh->vertices()[thirdVertexID];
    PRECICE_CHECK(utils::unique_elements(utils::make_array(vertices[0]->getCoords(),
                                                           vertices[1]->getCoords(), vertices[2]->getCoords())),
                  "setMeshTriangleWithEdges() was called with vertices located at identical coordinates (IDs: {}, {}, {}).",
                  firstVertexID, secondVertexID, thirdVertexID);
    mesh::Edge *edges[3];
    edges[0] = &mesh->createUniqueEdge(*vertices[0], *vertices[1]);
    edges[1] = &mesh->createUniqueEdge(*vertices[1], *vertices[2]);
    edges[2] = &mesh->createUniqueEdge(*vertices[2], *vertices[0]);

    mesh->createTriangle(*edges[0], *edges[1], *edges[2]);
  }
}

void SolverInterfaceImpl::setMeshQuad(
    MeshID meshID,
    int    firstEdgeID,
    int    secondEdgeID,
    int    thirdEdgeID,
    int    fourthEdgeID)
{
  PRECICE_TRACE(meshID, firstEdgeID, secondEdgeID, thirdEdgeID,
                fourthEdgeID);
  PRECICE_CHECK(_dimensions == 3, "setMeshQuad is only possible for 3D cases. "
                                  "Please set the dimension to 3 in the preCICE configuration file.");
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidEdgeID;
    PRECICE_CHECK(mesh->isValidEdgeID(firstEdgeID), errorInvalidEdgeID(firstEdgeID));
    PRECICE_CHECK(mesh->isValidEdgeID(secondEdgeID), errorInvalidEdgeID(secondEdgeID));
    PRECICE_CHECK(mesh->isValidEdgeID(thirdEdgeID), errorInvalidEdgeID(thirdEdgeID));
    PRECICE_CHECK(mesh->isValidEdgeID(fourthEdgeID), errorInvalidEdgeID(fourthEdgeID));

    PRECICE_CHECK(utils::unique_elements(utils::make_array(firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID)),
                  "The four edge ID's are not unique. Please check that the edges that form the quad are correct.");

    auto chain = mesh::asChain(utils::make_array(
        &mesh->edges()[firstEdgeID], &mesh->edges()[secondEdgeID],
        &mesh->edges()[thirdEdgeID], &mesh->edges()[fourthEdgeID]));
    PRECICE_CHECK(chain.connected, "The four edges are not connect. Please check that the edges that form the quad are correct.");

    auto coords = mesh::coordsFor(chain.vertices);
    PRECICE_CHECK(utils::unique_elements(coords),
                  "The four vertices that form the quad are not unique. "
                  "The resulting shape may be a point, line or triangle."
                  "Please check that the adapter sends the four unique vertices that form the quad, or that the mesh on the interface "
                  "is composed of planar quads.");

    auto convexity = math::geometry::isConvexQuad(coords);
    PRECICE_CHECK(convexity.convex,
                  "The given quad is not convex. "
                  "Please check that the adapter send the four correct vertices or that the interface is composed of planar quads.");

    // Use the shortest diagonal to split the quad into 2 triangles.
    // The diagonal to be used with edges (1, 2) and (0, 3) of the chain
    double distance1 = (coords[0] - coords[2]).norm();
    // The diagonal to be used with edges (0, 1) and (2, 3) of the chain
    double distance2 = (coords[1] - coords[3]).norm();

    // The new edge, e[4], is the shortest diagonal of the quad
    if (distance1 <= distance2) {
      auto &diag = mesh->createUniqueEdge(*chain.vertices[0], *chain.vertices[2]);
      mesh->createTriangle(*chain.edges[3], *chain.edges[0], diag);
      mesh->createTriangle(*chain.edges[1], *chain.edges[2], diag);
    } else {
      auto &diag = mesh->createUniqueEdge(*chain.vertices[1], *chain.vertices[3]);
      mesh->createTriangle(*chain.edges[0], *chain.edges[1], diag);
      mesh->createTriangle(*chain.edges[2], *chain.edges[3], diag);
    }
  }
}

void SolverInterfaceImpl::setMeshQuadWithEdges(
    MeshID meshID,
    int    firstVertexID,
    int    secondVertexID,
    int    thirdVertexID,
    int    fourthVertexID)
{
  PRECICE_TRACE(meshID, firstVertexID,
                secondVertexID, thirdVertexID, fourthVertexID);
  PRECICE_CHECK(_dimensions == 3, "setMeshQuadWithEdges is only possible for 3D cases."
                                  " Please set the dimension to 3 in the preCICE configuration file.");
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext &context = _accessor->usedMeshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    PRECICE_ASSERT(context.mesh);
    mesh::Mesh &mesh = *(context.mesh);
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh.isValidVertexID(firstVertexID), errorInvalidVertexID(firstVertexID));
    PRECICE_CHECK(mesh.isValidVertexID(secondVertexID), errorInvalidVertexID(secondVertexID));
    PRECICE_CHECK(mesh.isValidVertexID(thirdVertexID), errorInvalidVertexID(thirdVertexID));
    PRECICE_CHECK(mesh.isValidVertexID(fourthVertexID), errorInvalidVertexID(fourthVertexID));

    auto vertexIDs = utils::make_array(firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
    PRECICE_CHECK(utils::unique_elements(vertexIDs), "The four vertex ID's are not unique. Please check that the vertices that form the quad are correct.");

    auto coords = mesh::coordsFor(mesh, vertexIDs);
    PRECICE_CHECK(utils::unique_elements(coords),
                  "The four vertices that form the quad are not unique. The resulting shape may be a point, line or triangle."
                  "Please check that the adapter sends the four unique vertices that form the quad, or that the mesh on the interface "
                  "is composed of quads. A mix of triangles and quads are not supported.");

    auto convexity = math::geometry::isConvexQuad(coords);
    PRECICE_CHECK(convexity.convex, "The given quad is not convex. "
                                    "Please check that the adapter send the four correct vertices or that the interface is composed of quads. "
                                    "A mix of triangles and quads are not supported.");
    auto reordered = utils::reorder_array(convexity.vertexOrder, mesh::vertexPtrsFor(mesh, vertexIDs));

    // Vertices are now in the order: V0-V1-V2-V3-V0.
    // The order now identifies all outer edges of the quad.
    auto &edge0 = mesh.createUniqueEdge(*reordered[0], *reordered[1]);
    auto &edge1 = mesh.createUniqueEdge(*reordered[1], *reordered[2]);
    auto &edge2 = mesh.createUniqueEdge(*reordered[2], *reordered[3]);
    auto &edge3 = mesh.createUniqueEdge(*reordered[3], *reordered[0]);

    // Use the shortest diagonal to split the quad into 2 triangles.
    // Vertices are now in V0-V1-V2-V3-V0 order. The new edge, e[4] is either 0-2 or 1-3
    double distance1 = (reordered[0]->getCoords() - reordered[2]->getCoords()).norm();
    double distance2 = (reordered[1]->getCoords() - reordered[3]->getCoords()).norm();

    // The new edge, e[4], is the shortest diagonal of the quad
    if (distance1 <= distance2) {
      auto &diag = mesh.createUniqueEdge(*reordered[0], *reordered[2]);
      mesh.createTriangle(edge0, edge1, diag);
      mesh.createTriangle(edge2, edge3, diag);
    } else {
      auto &diag = mesh.createUniqueEdge(*reordered[1], *reordered[3]);
      mesh.createTriangle(edge3, edge0, diag);
      mesh.createTriangle(edge1, edge2, diag);
    }
  }
}

void SolverInterfaceImpl::mapWriteDataFrom(
    int fromMeshID)
{
  PRECICE_TRACE(fromMeshID);
  PRECICE_VALIDATE_MESH_ID(fromMeshID);
  impl::MeshContext &context = _accessor->usedMeshContext(fromMeshID);

  PRECICE_CHECK(not context.fromMappingContexts.empty(),
                "You attempt to \"mapWriteDataFrom\" mesh {}, but there is no mapping from this mesh configured. Maybe you don't want to call this function at all or you forgot to configure the mapping.",
                context.mesh->getName());

  double time = _couplingScheme->getTime();
  performDataActions({action::Action::WRITE_MAPPING_PRIOR}, time, 0, 0, 0);

  for (impl::MappingContext &mappingContext : context.fromMappingContexts) {
    if (not mappingContext.mapping->hasComputedMapping()) {
      PRECICE_DEBUG("Compute mapping from mesh \"{}\"", context.mesh->getName());
      mappingContext.mapping->computeMapping();
    }
    for (impl::DataContext &context : _accessor->writeDataContexts()) {

      if (context.mesh->getID() != fromMeshID) {
        continue;
      }

      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      PRECICE_DEBUG("Map data \"{}\" from mesh \"{}\"", context.fromData->getName(), context.mesh->getName());
      PRECICE_ASSERT(mappingContext.mapping == context.mappingContext.mapping);
      mappingContext.mapping->map(context.fromData->getID(), context.toData->getID());
    }
    mappingContext.hasMappedData = true;
  }
  performDataActions({action::Action::WRITE_MAPPING_POST}, time, 0, 0, 0);
}

void SolverInterfaceImpl::mapReadDataTo(
    int toMeshID)
{
  PRECICE_TRACE(toMeshID);
  PRECICE_VALIDATE_MESH_ID(toMeshID);
  impl::MeshContext &context = _accessor->usedMeshContext(toMeshID);

  PRECICE_CHECK(not context.toMappingContexts.empty(),
                "You attempt to \"mapReadDataTo\" mesh {}, but there is no mapping to this mesh configured. Maybe you don't want to call this function at all or you forgot to configure the mapping.",
                context.mesh->getName());

  double time = _couplingScheme->getTime();
  performDataActions({action::Action::READ_MAPPING_PRIOR}, time, 0, 0, 0);

  for (impl::MappingContext &mappingContext : context.toMappingContexts) {
    if (not mappingContext.mapping->hasComputedMapping()) {
      PRECICE_DEBUG("Compute mapping from mesh \"{}\"", context.mesh->getName());
      mappingContext.mapping->computeMapping();
    }
    for (impl::DataContext &context : _accessor->readDataContexts()) {
      if (context.mesh->getID() != toMeshID) {
        continue;
      }
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      PRECICE_DEBUG("Map data \"{}\" to mesh \"{}\"", context.fromData->getName(), context.mesh->getName());
      PRECICE_ASSERT(mappingContext.mapping == context.mappingContext.mapping);
      mappingContext.mapping->map(context.fromData->getID(), context.toData->getID());
      PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, context.toData->values()));
    }
    mappingContext.hasMappedData = true;
  }
  performDataActions({action::Action::READ_MAPPING_POST}, time, 0, 0, 0);
}

void SolverInterfaceImpl::writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_TRACE(dataID, size);
  PRECICE_CHECK(_state != State::Finalized, "writeBlockVectorData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_WRITE(dataID);
  if (size == 0)
    return;
  PRECICE_CHECK(valueIndices != nullptr, "writeBlockVectorData() was called with valueIndices == nullptr");
  PRECICE_CHECK(values != nullptr, "writeBlockVectorData() was called with values == nullptr");
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.fromData != nullptr);
  mesh::Data &data = *context.fromData;
  PRECICE_CHECK(data.getDimensions() == _dimensions,
                "You cannot call writeBlockVectorData on the scalar data type \"{0}\". Use writeBlockScalarData or change the data type for \"{0}\" to vector.",
                data.getName());
  PRECICE_VALIDATE_DATA(values, size * _dimensions);

  auto &     valuesInternal = data.values();
  const auto vertexCount    = valuesInternal.size() / data.getDimensions();
  for (int i = 0; i < size; i++) {
    const auto valueIndex = valueIndices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot write data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  data.getName(), valueIndex);
    const int offsetInternal = valueIndex * _dimensions;
    const int offset         = i * _dimensions;
    for (int dim = 0; dim < _dimensions; dim++) {
      PRECICE_ASSERT(offset + dim < valuesInternal.size(),
                     offset + dim, valuesInternal.size());
      valuesInternal[offsetInternal + dim] = values[offset + dim];
    }
  }
}

void SolverInterfaceImpl::writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *value)
{
  PRECICE_TRACE(dataID, valueIndex);
  PRECICE_CHECK(_state != State::Finalized, "writeVectorData(...) cannot be called before finalize().");
  PRECICE_REQUIRE_DATA_WRITE(dataID);
  PRECICE_DEBUG("value = {}", Eigen::Map<const Eigen::VectorXd>(value, _dimensions).format(utils::eigenio::debug()));
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.fromData != nullptr);
  mesh::Data &data = *context.fromData;
  PRECICE_CHECK(data.getDimensions() == _dimensions,
                "You cannot call writeVectorData on the scalar data type \"{0}\". Use writeScalarData or change the data type for \"{0}\" to vector.",
                data.getName());
  PRECICE_VALIDATE_DATA(value, _dimensions);

  auto &     values      = data.values();
  const auto vertexCount = values.size() / data.getDimensions();
  PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                "Cannot write data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                data.getName(), valueIndex);
  const int offset = valueIndex * _dimensions;
  for (int dim = 0; dim < _dimensions; dim++) {
    values[offset + dim] = value[dim];
  }
}

void SolverInterfaceImpl::writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_TRACE(dataID, size);
  PRECICE_CHECK(_state != State::Finalized, "writeBlockScalarData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_WRITE(dataID);
  if (size == 0)
    return;
  PRECICE_CHECK(valueIndices != nullptr, "writeBlockScalarData() was called with valueIndices == nullptr");
  PRECICE_CHECK(values != nullptr, "writeBlockScalarData() was called with values == nullptr");
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.fromData != nullptr);
  mesh::Data &data = *context.fromData;
  PRECICE_CHECK(data.getDimensions() == 1,
                "You cannot call writeBlockScalarData on the vector data type \"{}\". Use writeBlockVectorData or change the data type for \"{}\" to scalar.",
                data.getName(), data.getName());
  PRECICE_VALIDATE_DATA(values, size);

  auto &     valuesInternal = data.values();
  const auto vertexCount    = valuesInternal.size() / data.getDimensions();
  for (int i = 0; i < size; i++) {
    const auto valueIndex = valueIndices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot write data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  data.getName(), valueIndex);
    valuesInternal[valueIndex] = values[i];
  }
}

void SolverInterfaceImpl::writeScalarData(
    int    dataID,
    int    valueIndex,
    double value)
{
  PRECICE_TRACE(dataID, valueIndex, value);
  PRECICE_CHECK(_state != State::Finalized, "writeScalarData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_WRITE(dataID);
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.fromData != nullptr);
  mesh::Data &data = *context.fromData;
  PRECICE_CHECK(valueIndex >= -1,
                "Invalid value index ({}) when writing scalar data. Value index must be >= 0. "
                "Please check the value index for {}",
                valueIndex, data.getName());
  PRECICE_CHECK(data.getDimensions() == 1,
                "You cannot call writeScalarData on the vector data type \"{0}\". "
                "Use writeVectorData or change the data type for \"{0}\" to scalar.",
                data.getName());
  PRECICE_VALIDATE_DATA(static_cast<double *>(&value), 1);

  auto &     values      = data.values();
  const auto vertexCount = values.size() / data.getDimensions();
  PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                "Cannot write data \"{}\" to invalid Vertex ID ({}). "
                "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                data.getName(), valueIndex);
  values[valueIndex] = value;
}

void SolverInterfaceImpl::readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  PRECICE_TRACE(dataID, size);
  PRECICE_CHECK(_state != State::Finalized, "readBlockVectorData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_READ(dataID);
  if (size == 0)
    return;
  PRECICE_CHECK(valueIndices != nullptr, "readBlockVectorData() was called with valueIndices == nullptr");
  PRECICE_CHECK(values != nullptr, "readBlockVectorData() was called with values == nullptr");
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.toData != nullptr);
  mesh::Data &data = *context.toData;
  PRECICE_CHECK(data.getDimensions() == _dimensions,
                "You cannot call readBlockVectorData on the scalar data type \"{0}\". "
                "Use readBlockScalarData or change the data type for \"{0}\" to vector.",
                data.getName());
  auto &     valuesInternal = data.values();
  const auto vertexCount    = valuesInternal.size() / data.getDimensions();
  for (int i = 0; i < size; i++) {
    const auto valueIndex = valueIndices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot read data \"{}\" to invalid Vertex ID ({}). "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  data.getName(), valueIndex);
    int offsetInternal = valueIndex * _dimensions;
    int offset         = i * _dimensions;
    for (int dim = 0; dim < _dimensions; dim++) {
      values[offset + dim] = valuesInternal[offsetInternal + dim];
    }
  }
}

void SolverInterfaceImpl::readVectorData(
    int     dataID,
    int     valueIndex,
    double *value) const
{
  PRECICE_TRACE(dataID, valueIndex);
  PRECICE_CHECK(_state != State::Finalized, "readVectorData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_READ(dataID);
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.toData != nullptr);
  mesh::Data &data = *context.toData;
  PRECICE_CHECK(valueIndex >= -1,
                "Invalid value index ( {} ) when reading vector data. Value index must be >= 0. "
                "Please check the value index for {}",
                valueIndex, data.getName());
  PRECICE_CHECK(data.getDimensions() == _dimensions,
                "You cannot call readVectorData on the scalar data type \"{0}\". Use readScalarData or change the data type for \"{0}\" to vector.",
                data.getName());
  auto &     values      = data.values();
  const auto vertexCount = values.size() / data.getDimensions();
  PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                "Cannot read data \"{}\" to invalid Vertex ID ({}). "
                "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                data.getName(), valueIndex);
  int offset = valueIndex * _dimensions;
  for (int dim = 0; dim < _dimensions; dim++) {
    value[dim] = values[offset + dim];
  }
  PRECICE_DEBUG("read value = {}", Eigen::Map<const Eigen::VectorXd>(value, _dimensions).format(utils::eigenio::debug()));
}

void SolverInterfaceImpl::readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  PRECICE_TRACE(dataID, size);
  PRECICE_CHECK(_state != State::Finalized, "readBlockScalarData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_READ(dataID);
  if (size == 0)
    return;
  PRECICE_CHECK(valueIndices != nullptr, "readBlockScalarData() was called with valueIndices == nullptr");
  PRECICE_CHECK(values != nullptr, "readBlockScalarData() was called with values == nullptr");
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.toData != nullptr);
  mesh::Data &data = *context.toData;
  PRECICE_CHECK(data.getDimensions() == 1,
                "You cannot call readBlockScalarData on the vector data type \"{0}\". "
                "Use readBlockVectorData or change the data type for \"{0}\" to scalar.",
                data.getName());
  auto &     valuesInternal = data.values();
  const auto vertexCount    = valuesInternal.size();

  for (int i = 0; i < size; i++) {
    const auto valueIndex = valueIndices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot read data \"{}\" to invalid Vertex ID ({}). "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  data.getName(), valueIndex);
    values[i] = valuesInternal[valueIndex];
  }
}

void SolverInterfaceImpl::readScalarData(
    int     dataID,
    int     valueIndex,
    double &value) const
{
  PRECICE_TRACE(dataID, valueIndex, value);
  PRECICE_CHECK(_state != State::Finalized, "readScalarData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_READ(dataID);
  DataContext &context = _accessor->dataContext(dataID);
  PRECICE_ASSERT(context.toData != nullptr);
  mesh::Data &data = *context.toData;
  PRECICE_CHECK(valueIndex >= -1,
                "Invalid value index ( {} ) when reading scalar data. Value index must be >= 0. "
                "Please check the value index for {}",
                valueIndex, data.getName());
  PRECICE_CHECK(data.getDimensions() == 1,
                "You cannot call readScalarData on the vector data type \"{}\". "
                "Use readVectorData or change the data type for \"{}\" to scalar.",
                data.getName());
  auto &     values      = data.values();
  const auto vertexCount = values.size();
  PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                "Cannot read data \"{}\" from invalid Vertex ID ({}). "
                "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                data.getName(), valueIndex);
  value = values[valueIndex];
  PRECICE_DEBUG("Read value = {}", value);
}

void SolverInterfaceImpl::setMeshAccessRegion(
    const int     meshID,
    const double *boundingBox) const
{
  PRECICE_TRACE(meshID);
  PRECICE_REQUIRE_MESH_USE(meshID);
  PRECICE_CHECK(_state != State::Finalized, "setMeshAccessRegion() cannot be called after finalize().")
  PRECICE_CHECK(_state != State::Initialized, "setMeshAccessRegion() needs to be called before initialize().");
  PRECICE_CHECK(!_accessRegionDefined, "setMeshAccessRegion may only be called once.");
  PRECICE_CHECK(boundingBox != nullptr, "setMeshAccessRegion was called with boundingBox == nullptr.");

  // Get the related mesh
  MeshContext & context = _accessor->meshContext(meshID);
  mesh::PtrMesh mesh(context.mesh);
  PRECICE_DEBUG("Define bounding box");
  // Transform bounds into a suitable format
  int                 dim = mesh->getDimensions();
  std::vector<double> bounds(dim * 2);

  for (int d = 0; d < dim; ++d) {
    // Check that min is lower or equal to max
    PRECICE_CHECK(boundingBox[2 * d] <= boundingBox[2 * d + 1], "Your bounding box is ill defined, i.e. it has a negative volume. The required format is [x_min, x_max...]");
    bounds[2 * d]     = boundingBox[2 * d];
    bounds[2 * d + 1] = boundingBox[2 * d + 1];
  }
  // Create a bounding box
  mesh::BoundingBox providedBoundingBox(bounds);
  // Expand the mesh associated bounding box
  mesh->expandBoundingBox(providedBoundingBox);
  // and set a flag so that we know the function was called
  _accessRegionDefined = true;
}

void SolverInterfaceImpl::getMeshVerticesAndIDs(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates) const
{
  PRECICE_TRACE(meshID, size);
  PRECICE_REQUIRE_MESH_USE(meshID);
  PRECICE_DEBUG("Get {} mesh vertices with IDs", size);
  if (size == 0)
    return;

  const MeshContext & context = _accessor->meshContext(meshID);
  const mesh::PtrMesh mesh(context.mesh);

  PRECICE_CHECK(ids != nullptr, "getMeshVerticesWithIDs() was called with ids == nullptr");
  PRECICE_CHECK(coordinates != nullptr, "getMeshVerticesWithIDs() was called with coordinates == nullptr");

  const auto &vertices = mesh->vertices();
  PRECICE_CHECK(size <= vertices.size(), "The queried size exceeds the number of available points.");

  Eigen::Map<Eigen::MatrixXd> posMatrix{
      coordinates, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};

  for (size_t i = 0; i < size; i++) {
    PRECICE_ASSERT(i < vertices.size(), i, vertices.size());
    ids[i]           = vertices[i].getID();
    posMatrix.col(i) = vertices[i].getCoords();
  }
}

void SolverInterfaceImpl::exportMesh(
    const std::string &filenameSuffix,
    int                exportType) const
{
  PRECICE_TRACE(filenameSuffix, exportType);
  // Export meshes
  //const ExportContext& context = _accessor->exportContext();
  for (const io::ExportContext &context : _accessor->exportContexts()) {
    PRECICE_DEBUG("Export type = {}", exportType);
    bool exportAll  = exportType == io::constants::exportAll();
    bool exportThis = context.exporter->getType() == exportType;
    if (exportAll || exportThis) {
      for (const MeshContext *meshContext : _accessor->usedMeshContexts()) {
        std::string name = meshContext->mesh->getName() + "-" + filenameSuffix;
        PRECICE_DEBUG("Exporting mesh to file \"{}\" at location \"{}\"", name, context.location);
        context.exporter->doExport(name, context.location, *(meshContext->mesh));
      }
    }
  }
}

void SolverInterfaceImpl::configureM2Ns(
    const m2n::M2NConfiguration::SharedPointer &config)
{
  PRECICE_TRACE();
  for (const auto &m2nTuple : config->m2ns()) {
    std::string comPartner("");
    bool        isRequesting = false;
    if (std::get<1>(m2nTuple) == _accessorName) {
      comPartner   = std::get<2>(m2nTuple);
      isRequesting = true;
    } else if (std::get<2>(m2nTuple) == _accessorName) {
      comPartner = std::get<1>(m2nTuple);
    }
    if (not comPartner.empty()) {
      for (const impl::PtrParticipant &participant : _participants) {
        if (participant->getName() == comPartner) {
          PRECICE_ASSERT(not utils::contained(comPartner, _m2ns), comPartner);
          PRECICE_ASSERT(std::get<0>(m2nTuple));

          _m2ns[comPartner] = [&] {
            m2n::BoundM2N bound;
            bound.m2n          = std::get<0>(m2nTuple);
            bound.localName    = _accessorName;
            bound.remoteName   = comPartner;
            bound.isRequesting = isRequesting;
            return bound;
          }();
        }
      }
    }
  }
}

void SolverInterfaceImpl::configurePartitions(
    const m2n::M2NConfiguration::SharedPointer &m2nConfig)
{
  PRECICE_TRACE();
  for (MeshContext *context : _accessor->usedMeshContexts()) {

    if (context->provideMesh) { // Accessor provides mesh
      PRECICE_CHECK(context->receiveMeshFrom.empty(),
                    "Participant \"{}\" cannot provide and receive mesh {}!",
                    _accessorName, context->mesh->getName());

      context->partition = partition::PtrPartition(new partition::ProvidedPartition(context->mesh));

      for (auto &receiver : _participants) {
        for (auto &receiverContext : receiver->usedMeshContexts()) {
          if (receiverContext->receiveMeshFrom == _accessorName && receiverContext->mesh->getName() == context->mesh->getName()) {
            // meshRequirement has to be copied from "from" to provide", since
            // mapping are only defined at "provide"
            if (receiverContext->meshRequirement > context->meshRequirement) {
              context->meshRequirement = receiverContext->meshRequirement;
            }

            m2n::PtrM2N m2n = m2nConfig->getM2N(receiver->getName(), _accessorName);
            m2n->createDistributedCommunication(context->mesh);
            context->partition->addM2N(m2n);
          }
        }
      }
      /// @todo support offset??

    } else { // Accessor receives mesh
      PRECICE_CHECK(not context->receiveMeshFrom.empty(),
                    "Participant \"{}\" must either provide or receive the mesh \"{}\". "
                    "Please define either a \"from\" or a \"provide\" attribute in the <use-mesh name=\"{}\"/> node of \"{}\".",
                    _accessorName, context->mesh->getName(), context->mesh->getName(), _accessorName);
      PRECICE_CHECK(not context->provideMesh,
                    "Participant \"{}\" cannot provide and receive mesh \"{}\" at the same time. "
                    "Please check your \"from\" and \"provide\" attributes in the <use-mesh name=\"{}\"/> node of \"{}\".",
                    _accessorName, context->mesh->getName(), context->mesh->getName(), _accessorName);
      std::string receiver(_accessorName);
      std::string provider(context->receiveMeshFrom);

      PRECICE_DEBUG("Receiving mesh from {}", provider);

      context->partition = partition::PtrPartition(new partition::ReceivedPartition(context->mesh, context->geoFilter, context->safetyFactor, context->allowDirectAccess));

      m2n::PtrM2N m2n = m2nConfig->getM2N(receiver, provider);
      m2n->createDistributedCommunication(context->mesh);
      context->partition->addM2N(m2n);
      for (const MappingContext &mappingContext : context->fromMappingContexts) {
        context->partition->addFromMapping(mappingContext.mapping);
      }
      for (const MappingContext &mappingContext : context->toMappingContexts) {
        context->partition->addToMapping(mappingContext.mapping);
      }
    }
  }
}

void SolverInterfaceImpl::compareBoundingBoxes()
{
  // sort meshContexts by name, for communication in right order.
  std::sort(_accessor->usedMeshContexts().begin(), _accessor->usedMeshContexts().end(),
            [](MeshContext const *const lhs, MeshContext const *const rhs) -> bool {
              return lhs->mesh->getName() < rhs->mesh->getName();
            });

  for (MeshContext *meshContext : _accessor->usedMeshContexts()) {
    if (meshContext->provideMesh) // provided meshes need their bounding boxes already for the re-partitioning
      meshContext->mesh->computeBoundingBox();
  }

  for (MeshContext *meshContext : _accessor->usedMeshContexts()) {
    meshContext->partition->compareBoundingBoxes();
  }
}

void SolverInterfaceImpl::computePartitions()
{
  //We need to do this in two loops: First, communicate the mesh and later compute the partition.
  //Originally, this was done in one loop. This however gave deadlock if two meshes needed to be communicated cross-wise.
  //Both loops need a different sorting

  auto &contexts = _accessor->usedMeshContexts();

  std::sort(contexts.begin(), contexts.end(),
            [](MeshContext const *const lhs, MeshContext const *const rhs) -> bool {
              return lhs->mesh->getName() < rhs->mesh->getName();
            });

  for (MeshContext *meshContext : contexts) {
    meshContext->partition->communicate();
  }

  // for two-level initialization, there is also still communication in partition::compute()
  // therefore, we cannot resort here.
  // @todo this hacky solution should be removed as part of #633
  bool resort = true;
  for (auto &m2nPair : _m2ns) {
    if (m2nPair.second.m2n->usesTwoLevelInitialization()) {
      resort = false;
      break;
    }
  }

  if (resort) {
    // pull provided meshes up front, to have them ready for the decomposition of the received meshes (for the mappings)
    std::stable_partition(contexts.begin(), contexts.end(),
                          [](MeshContext const *const meshContext) -> bool {
                            return meshContext->provideMesh;
                          });
  }

  for (MeshContext *meshContext : contexts) {
    meshContext->partition->compute();
    if (not meshContext->provideMesh) { // received mesh can only compute their bounding boxes here
      meshContext->mesh->computeBoundingBox();
    }
    meshContext->mesh->allocateDataValues();
  }
}

void SolverInterfaceImpl::computeMappings(utils::ptr_vector<MappingContext> contexts, const std::string &mappingType)
{
  PRECICE_TRACE();
  using namespace mapping;
  MappingConfiguration::Timing timing;
  for (impl::MappingContext &context : contexts) {
    timing      = context.timing;
    bool mapNow = timing == MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == MappingConfiguration::INITIAL;
    bool hasComputed = context.mapping->hasComputedMapping();
    if (mapNow && not hasComputed) {
      PRECICE_INFO("Compute \"{}\" mapping from mesh \"{}\" to mesh \"{}\".",
                   mappingType, _accessor->meshContext(context.fromMeshID).mesh->getName(), _accessor->meshContext(context.toMeshID).mesh->getName());

      context.mapping->computeMapping();
    }
  }
}

void SolverInterfaceImpl::mapData(utils::ptr_vector<DataContext> contexts, const std::string &mappingType)
{
  PRECICE_TRACE();
  using namespace mapping;
  MappingConfiguration::Timing timing;
  for (impl::DataContext &context : contexts) {
    timing          = context.mappingContext.timing;
    bool hasMapping = context.mappingContext.mapping.get() != nullptr;
    bool hasMapped  = context.mappingContext.hasMappedData;
    bool mapNow     = timing == MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == MappingConfiguration::INITIAL;
    if (mapNow && hasMapping && (not hasMapped)) {
      int inDataID  = context.fromData->getID();
      int outDataID = context.toData->getID();
      PRECICE_DEBUG("Map \"{}\" data \"{}\" from mesh \"{}\"",
                    mappingType, context.fromData->getName(), context.mesh->getName());
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      PRECICE_DEBUG("Map from dataID {} to dataID: {}", inDataID, outDataID);
      context.mappingContext.mapping->map(inDataID, outDataID);
      PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, context.toData->values()));
    }
  }
}

void SolverInterfaceImpl::clearMappings(utils::ptr_vector<MappingContext> contexts)
{
  PRECICE_TRACE();
  // Clear non-stationary, non-incremental mappings
  using namespace mapping;
  for (impl::MappingContext &context : contexts) {
    bool isStationary = context.timing == MappingConfiguration::INITIAL;
    if (not isStationary) {
      context.mapping->clear();
    }
    context.hasMappedData = false;
  }
}

void SolverInterfaceImpl::mapWrittenData()
{
  PRECICE_TRACE();
  computeMappings(_accessor->writeMappingContexts(), "write");
  mapData(_accessor->writeDataContexts(), "write");
  clearMappings(_accessor->writeMappingContexts());
}

void SolverInterfaceImpl::mapReadData()
{
  PRECICE_TRACE();
  computeMappings(_accessor->readMappingContexts(), "read");
  mapData(_accessor->readDataContexts(), "read");
  clearMappings(_accessor->readMappingContexts());
}

void SolverInterfaceImpl::performDataActions(
    const std::set<action::Action::Timing> &timings,
    double                                  time,
    double                                  timeStepSize,
    double                                  computedTimeWindowPart,
    double                                  timeWindowSize)
{
  PRECICE_TRACE();
  for (action::PtrAction &action : _accessor->actions()) {
    if (timings.find(action->getTiming()) != timings.end()) {
      action->performAction(time, timeStepSize, computedTimeWindowPart, timeWindowSize);
    }
  }
}

void SolverInterfaceImpl::handleExports()
{
  PRECICE_TRACE();
  //timesteps was already incremented before
  int timesteps = _couplingScheme->getTimeWindows() - 1;

  for (const io::ExportContext &context : _accessor->exportContexts()) {
    if (_couplingScheme->isTimeWindowComplete() || context.everyIteration) {
      if (context.everyNTimeWindows != -1) {
        if (timesteps % context.everyNTimeWindows == 0) {
          if (context.everyIteration) {
            std::ostringstream everySuffix;
            everySuffix << _accessorName << ".it" << _numberAdvanceCalls;
            exportMesh(everySuffix.str());
          }
          std::ostringstream suffix;
          suffix << _accessorName << ".dt" << _couplingScheme->getTimeWindows() - 1;
          exportMesh(suffix.str());
        }
      }
    }
  }

  if (_couplingScheme->isTimeWindowComplete()) {
    // Export watch point data
    for (const PtrWatchPoint &watchPoint : _accessor->watchPoints()) {
      watchPoint->exportPointData(_couplingScheme->getTime());
    }
    for (const PtrWatchIntegral &watchIntegral : _accessor->watchIntegrals()) {
      watchIntegral->exportIntegralData(_couplingScheme->getTime());
    }
  }
}

void SolverInterfaceImpl::resetWrittenData()
{
  PRECICE_TRACE();
  for (DataContext &context : _accessor->writeDataContexts()) {
    context.fromData->toZero();
    if (context.toData != context.fromData) {
      context.toData->toZero();
    }
  }
}

PtrParticipant SolverInterfaceImpl::determineAccessingParticipant(
    const config::SolverInterfaceConfiguration &config)
{
  const auto &partConfig = config.getParticipantConfiguration();
  for (const PtrParticipant &participant : partConfig->getParticipants()) {
    if (participant->getName() == _accessorName) {
      return participant;
    }
  }
  PRECICE_ERROR("This participant's name, which was specified in the constructor of the preCICE interface as \"{}\", "
                "is not defined in the preCICE configuration. "
                "Please double-check the correct spelling.",
                _accessorName);
}

void SolverInterfaceImpl::initializeMasterSlaveCommunication()
{
  PRECICE_TRACE();

  Event e("com.initializeMasterSlaveCom", precice::syncMode);
  utils::MasterSlave::_communication->connectMasterSlaves(
      _accessorName, "MasterSlaves",
      _accessorProcessRank, _accessorCommunicatorSize);
}

void SolverInterfaceImpl::syncTimestep(double computedTimestepLength)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel());
  if (utils::MasterSlave::isSlave()) {
    utils::MasterSlave::_communication->send(computedTimestepLength, 0);
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isMaster());
    for (Rank rankSlave : utils::MasterSlave::allSlaves()) {
      double dt;
      utils::MasterSlave::_communication->receive(dt, rankSlave);
      PRECICE_CHECK(math::equals(dt, computedTimestepLength),
                    "Found ambiguous values for the timestep length passed to preCICE in \"advance\". On rank {}, the value is {}, while on rank 0, the value is {}.",
                    rankSlave, dt, computedTimestepLength);
    }
  }
}

void SolverInterfaceImpl::closeCommunicationChannels(CloseChannels close)
{
  // Apply some final ping-pong to synch solver that run e.g. with a uni-directional coupling only
  // afterwards close connections
  PRECICE_INFO("Synchronize participants and close {}communication channels",
               (close == CloseChannels::Distributed ? "distributed " : ""));
  std::string ping = "ping";
  std::string pong = "pong";
  for (auto &iter : _m2ns) {
    auto bm2n = iter.second;
    if (not utils::MasterSlave::isSlave()) {
      PRECICE_DEBUG("Synchronizing Master with {}", bm2n.remoteName);
      if (bm2n.isRequesting) {
        bm2n.m2n->getMasterCommunication()->send(ping, 0);
        std::string receive = "init";
        bm2n.m2n->getMasterCommunication()->receive(receive, 0);
        PRECICE_ASSERT(receive == pong);
      } else {
        std::string receive = "init";
        bm2n.m2n->getMasterCommunication()->receive(receive, 0);
        PRECICE_ASSERT(receive == ping);
        bm2n.m2n->getMasterCommunication()->send(pong, 0);
      }
    }
    if (close == CloseChannels::Distributed) {
      PRECICE_DEBUG("Closing distributed communication with {}", bm2n.remoteName);
      bm2n.m2n->closeDistributedConnections();
    } else {
      PRECICE_DEBUG("Closing communication with {}", bm2n.remoteName);
      bm2n.m2n->closeConnection();
    }
  }
}

const mesh::Mesh &SolverInterfaceImpl::mesh(const std::string &meshName) const
{
  PRECICE_TRACE(meshName);
  return *_accessor->usedMeshContext(meshName).mesh;
}

} // namespace impl
} // namespace precice

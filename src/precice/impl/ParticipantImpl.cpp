#include <Eigen/Core>
#include <Eigen/src/Core/util/Meta.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <filesystem>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <tuple>
#include <utility>

#include "ParticipantImpl.hpp"
#include "action/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
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
#include "mapping/device/Ginkgo.hpp"
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
#include "precice/impl/CommonErrorMessages.hpp"
#include "precice/impl/MappingContext.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/ReadDataContext.hpp"
#include "precice/impl/Types.hpp"
#include "precice/impl/ValidationMacros.hpp"
#include "precice/impl/WatchIntegral.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "precice/impl/WriteDataContext.hpp"
#include "precice/impl/versions.hpp"
#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"
#include "profiling/config/ProfilingConfiguration.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/EigenIO.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/algorithm.hpp"
#include "utils/assertion.hpp"
#include "xml/XMLTag.hpp"

using precice::profiling::Event;

namespace precice::impl {

ParticipantImpl::ParticipantImpl(
    std::string_view      participantName,
    std::string_view      configurationFileName,
    int                   solverProcessIndex,
    int                   solverProcessSize,
    std::optional<void *> communicator)
    : _accessorName(participantName),
      _accessorProcessRank(solverProcessIndex),
      _accessorCommunicatorSize(solverProcessSize)
{

  PRECICE_CHECK(!_accessorName.empty(),
                "This participant's name is an empty string. "
                "When constructing a preCICE interface you need to pass the name of the "
                "participant as first argument to the constructor.");
  logging::setParticipant(_accessorName);

  PRECICE_CHECK(!communicator || communicator.value() != nullptr,
                "Passing \"nullptr\" as \"communicator\" to Participant constructor is not allowed. "
                "Please use the Participant constructor without the \"communicator\" argument, if you don't want to pass an MPI communicator.");
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

  profiling::EventRegistry::instance().initialize(_accessorName, _accessorProcessRank, _accessorCommunicatorSize);
  profiling::applyDefaults();
  Event e("construction", profiling::Fundamental);

  // Set the global communicator to the passed communicator.
  // This is a noop if preCICE is not configured with MPI.
#ifndef PRECICE_NO_MPI
  Event e3("com.initializeMPI", profiling::Fundamental);
  if (communicator.has_value()) {
    auto commptr = static_cast<utils::Parallel::Communicator *>(communicator.value());
    utils::Parallel::initializeOrDetectMPI(*commptr);
  } else {
    utils::Parallel::initializeOrDetectMPI();
  }

  {
    const auto currentRank = utils::Parallel::current()->rank();
    PRECICE_CHECK(_accessorProcessRank == currentRank,
                  "The solver process index given in the preCICE interface constructor({}) does not match the rank of the passed MPI communicator ({}).",
                  _accessorProcessRank, currentRank);
    const auto currentSize = utils::Parallel::current()->size();
    PRECICE_CHECK(_accessorCommunicatorSize == currentSize,
                  "The solver process size given in the preCICE interface constructor({}) does not match the size of the passed MPI communicator ({}).",
                  _accessorCommunicatorSize, currentSize);
  }
  e3.stop();
#endif

  Event e1("configure", profiling::Fundamental);
  configure(configurationFileName);
  e1.stop();

  // Backend settings have been configured
  Event e2("startProfilingBackend");
  profiling::EventRegistry::instance().startBackend();
  e2.stop();

  PRECICE_DEBUG("Initialize intra-participant communication");
  if (utils::IntraComm::isParallel()) {
    initializeIntraCommunication();
  }

  e.stop();
  _solverInitEvent = std::make_unique<profiling::Event>("solver.initialize", profiling::Fundamental, profiling::Synchronize);
}

ParticipantImpl::~ParticipantImpl()
{
  if (_state != State::Finalized) {
    PRECICE_INFO("Implicitly finalizing in destructor");
    finalize();
  }
}

void ParticipantImpl::configure(
    std::string_view configurationFileName)
{

  config::Configuration config;
  logging::setMPIRank(utils::Parallel::current()->rank());
  xml::ConfigurationContext context{
      _accessorName,
      _accessorProcessRank,
      _accessorCommunicatorSize};
  xml::configure(config.getXMLTag(), context, configurationFileName);
  if (_accessorProcessRank == 0) {
    PRECICE_INFO("This is preCICE version {}", PRECICE_VERSION);
    PRECICE_INFO("Revision info: {}", precice::preciceRevision);
    PRECICE_INFO("Build type: "
#ifndef NDEBUG
                 "Debug"
#else // NDEBUG
                 "Release"
#ifndef PRECICE_NO_DEBUG_LOG
                 " + debug log"
#else
                 " (without debug log)"
#endif
#ifndef PRECICE_NO_TRACE_LOG
                 " + trace log"
#endif
#ifndef PRECICE_NO_ASSERTIONS
                 " + assertions"
#endif
#endif // NDEBUG
    );
    try {
      PRECICE_INFO("Working directory \"{}\"", std::filesystem::current_path().string());
    } catch (std::filesystem::filesystem_error &fse) {
      PRECICE_INFO("Working directory unknown due to error \"{}\"", fse.what());
    }
    PRECICE_INFO("Configuring preCICE with configuration \"{}\"", configurationFileName);
    PRECICE_INFO("I am participant \"{}\"", _accessorName);
  }

  PRECICE_TRACE();

  _meshLock.clear();

  _allowsExperimental = config.allowsExperimental();
  _allowsRemeshing    = config.allowsRemeshing();
  _waitInFinalize     = config.waitInFinalize();
  _accessor           = determineAccessingParticipant(config);
  _accessor->setMeshIdManager(config.getMeshConfiguration()->extractMeshIdManager());

  PRECICE_ASSERT(_accessorCommunicatorSize == 1 || _accessor->useIntraComm(),
                 "A parallel participant needs an intra-participant communication");
  PRECICE_CHECK(not(_accessorCommunicatorSize == 1 && _accessor->useIntraComm()),
                "You cannot use an intra-participant communication with a serial participant. "
                "If you do not know exactly what an intra-participant communication is and why you want to use it "
                "you probably just want to remove the intraComm tag from the preCICE configuration.");

  utils::IntraComm::configure(_accessorProcessRank, _accessorCommunicatorSize);

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
    _meshLock.add(meshContext->mesh->getName(), false);
  }
}

void ParticipantImpl::initialize()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Finalized, "initialize() cannot be called after finalize().");
  PRECICE_CHECK(_state != State::Initialized, "initialize() may only be called once.");
  PRECICE_ASSERT(not _couplingScheme->isInitialized());

  bool failedToInitialize = _couplingScheme->isActionRequired(cplscheme::CouplingScheme::Action::InitializeData) && not _couplingScheme->isActionFulfilled(cplscheme::CouplingScheme::Action::InitializeData);
  PRECICE_CHECK(not failedToInitialize,
                "Initial data has to be written to preCICE before calling initialize(). "
                "After defining your mesh, call requiresInitialData() to check if the participant is required to write initial data using the writeData() function.");

  _solverInitEvent.reset();
  Event e("initialize", profiling::Fundamental, profiling::Synchronize);

  setupCommunication();
  setupWatcher();

  _meshLock.lockAll();

  for (auto &context : _accessor->writeDataContexts()) {
    const double startTime = 0.0;
    context.storeBufferedData(startTime);
  }

  mapInitialWrittenData();
  performDataActions({action::Action::WRITE_MAPPING_POST});

  PRECICE_DEBUG("Initialize coupling schemes");
  Event e1("initalizeCouplingScheme", profiling::Fundamental);
  _couplingScheme->initialize();
  e1.stop();

  mapInitialReadData();
  performDataActions({action::Action::READ_MAPPING_POST});

  handleExports(ExportTiming::Initial);

  resetWrittenData();

  e.stop();

  _state = State::Initialized;
  PRECICE_INFO(_couplingScheme->printCouplingState());
  _solverAdvanceEvent = std::make_unique<profiling::Event>("solver.advance", profiling::Fundamental, profiling::Synchronize);
}

void ParticipantImpl::reinitialize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_allowsRemeshing);
  PRECICE_INFO("Reinitializing Participant");
  Event e("reinitialize", profiling::Fundamental);
  closeCommunicationChannels(CloseChannels::Distributed);

  setupCommunication();
  setupWatcher();

  PRECICE_DEBUG("Reinitialize coupling schemes");
  _couplingScheme->reinitialize();
}

void ParticipantImpl::setupCommunication()
{
  PRECICE_TRACE();

  // TODO only preprocess changed meshes
  PRECICE_DEBUG("Preprocessing provided meshes");
  for (MeshContext *meshContext : _accessor->usedMeshContexts()) {
    if (meshContext->provideMesh) {
      auto &mesh = *(meshContext->mesh);
      Event e("preprocess." + mesh.getName());
      meshContext->mesh->preprocess();
    }
  }

  // Setup communication

  PRECICE_INFO("Setting up primary communication to coupling partner/s");
  Event e2("connectPrimaries");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n       = m2nPair.second;
    bool  requesting = bm2n.isRequesting;
    if (bm2n.m2n->isConnected()) {
      PRECICE_DEBUG("Primary connection {} {} already connected.", (requesting ? "from" : "to"), bm2n.remoteName);
    } else {
      PRECICE_DEBUG((requesting ? "Awaiting primary connection from {}" : "Establishing primary connection to {}"), bm2n.remoteName);
      bm2n.prepareEstablishment();
      bm2n.connectPrimaryRanks();
      PRECICE_DEBUG("Established primary connection {} {}", (requesting ? "from " : "to "), bm2n.remoteName);
    }
  }
  e2.stop();

  PRECICE_INFO("Primary ranks are connected");

  Event e3("repartitioning");
  compareBoundingBoxes();

  PRECICE_INFO("Setting up preliminary secondary communication to coupling partner/s");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n = m2nPair.second;
    bm2n.preConnectSecondaryRanks();
  }

  computePartitions();
  e3.stop();

  PRECICE_INFO("Setting up secondary communication to coupling partner/s");
  Event e4("connectSecondaries");
  for (auto &m2nPair : _m2ns) {
    auto &bm2n = m2nPair.second;
    bm2n.connectSecondaryRanks();
    PRECICE_DEBUG("Established secondary connection {} {}", (bm2n.isRequesting ? "from " : "to "), bm2n.remoteName);
  }
  PRECICE_INFO("Secondary ranks are connected");

  for (auto &m2nPair : _m2ns) {
    m2nPair.second.cleanupEstablishment();
  }
}

void ParticipantImpl::setupWatcher()
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Initialize watchpoints");
  for (PtrWatchPoint &watchPoint : _accessor->watchPoints()) {
    watchPoint->initialize();
  }
  for (PtrWatchIntegral &watchIntegral : _accessor->watchIntegrals()) {
    watchIntegral->initialize();
  }
}

void ParticipantImpl::advance(
    double computedTimeStepSize)
{

  PRECICE_TRACE(computedTimeStepSize);

  // Events for the solver time, stopped when we enter, restarted when we leave advance
  PRECICE_ASSERT(_solverAdvanceEvent, "The advance event is created in initialize");
  _solverAdvanceEvent->stop();

  Event e("advance", profiling::Fundamental, profiling::Synchronize);

  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before advance().");
  PRECICE_CHECK(_state != State::Finalized, "advance() cannot be called after finalize().");
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before advance().");
  PRECICE_ASSERT(_couplingScheme->isInitialized());
  PRECICE_CHECK(isCouplingOngoing(), "advance() cannot be called when isCouplingOngoing() returns false.");
  PRECICE_CHECK(!math::equals(computedTimeStepSize, 0.0), "advance() cannot be called with a time step size of 0.");
  PRECICE_CHECK(computedTimeStepSize > 0.0, "advance() cannot be called with a negative time step size {}.", computedTimeStepSize);
  _numberAdvanceCalls++;

#ifndef NDEBUG
  PRECICE_DEBUG("Synchronize time step size");
  if (utils::IntraComm::isParallel()) {
    syncTimestep(computedTimeStepSize);
  }
#endif

  // Update the coupling scheme time state. Necessary to get correct remainder.
  const bool isAtWindowEnd = _couplingScheme->addComputedTime(computedTimeStepSize);

  if (_allowsRemeshing) {
    if (isAtWindowEnd) {
      int totalMeshChanges = getTotalMeshChanges();
      if (reinitHandshake(totalMeshChanges)) {
        reinitialize();
      }
    } else {
      PRECICE_CHECK(_meshLock.checkAll(), "The time window needs to end after remeshing.");
    }
  }

  const double timeSteppedTo = _couplingScheme->getTime();
  const auto   dataToReceive = _couplingScheme->implicitDataToReceive();

  handleDataBeforeAdvance(isAtWindowEnd, timeSteppedTo);

  advanceCouplingScheme();

  // In clase if an implicit scheme, this may be before timeSteppedTo
  const double timeAfterAdvance   = _couplingScheme->getTime();
  const bool   timeWindowComplete = _couplingScheme->isTimeWindowComplete();

  handleDataAfterAdvance(isAtWindowEnd, timeWindowComplete, timeSteppedTo, timeAfterAdvance, dataToReceive);

  PRECICE_INFO(_couplingScheme->printCouplingState());

  PRECICE_DEBUG("Mapped {} samples in write mappings and {} samples in read mappings",
                _executedWriteMappings, _executedReadMappings);

  _meshLock.lockAll();

  e.stop();
  _solverAdvanceEvent->start();
}

void ParticipantImpl::handleDataBeforeAdvance(bool reachedTimeWindowEnd, double timeSteppedTo)
{
  if (reachedTimeWindowEnd || _couplingScheme->requiresSubsteps()) {
    samplizeWriteData(timeSteppedTo);
  }
  resetWrittenData();

  // Reset mapping counters here to cover subcycling
  _executedReadMappings  = 0;
  _executedWriteMappings = 0;

  if (reachedTimeWindowEnd) {
    mapWrittenData(_couplingScheme->getTimeWindowStart());
    performDataActions({action::Action::WRITE_MAPPING_POST});
  }
}

void ParticipantImpl::handleDataAfterAdvance(bool reachedTimeWindowEnd, bool isTimeWindowComplete, double timeSteppedTo, double timeAfterAdvance, const cplscheme::ImplicitData &receivedData)
{
  if (!reachedTimeWindowEnd) {
    // We are subcycling
    return;
  }

  if (isTimeWindowComplete) {
    // Move to next time window
    PRECICE_ASSERT(math::greaterEquals(timeAfterAdvance, timeSteppedTo), "We must have stayed or moved forwards in time (min-time-step-size).", timeAfterAdvance, timeSteppedTo);

    // As we move forward, there may now be old samples lying around
    trimOldDataBefore(_couplingScheme->getTimeWindowStart());
  } else {
    // We are iterating
    PRECICE_ASSERT(math::greater(timeSteppedTo, timeAfterAdvance), "We must have moved back in time!");

    trimSendDataAfter(timeAfterAdvance);
  }

  if (reachedTimeWindowEnd) {
    trimReadMappedData(timeAfterAdvance, isTimeWindowComplete, receivedData);
    mapReadData();
    performDataActions({action::Action::READ_MAPPING_POST});
  }

  // Required for implicit coupling
  for (auto &context : _accessor->readDataContexts()) {
    context.invalidateMappingCache();
  }

  for (auto &context : _accessor->writeDataContexts()) {
    context.invalidateMappingCache();
  }

  if (isTimeWindowComplete) {
    // Reset initial guesses for iterative mappings
    for (auto &context : _accessor->readDataContexts()) {
      context.resetInitialGuesses();
    }
    for (auto &context : _accessor->writeDataContexts()) {
      context.resetInitialGuesses();
    }
  }

  handleExports(ExportTiming::Advance);
}

void ParticipantImpl::samplizeWriteData(double time)
{
  // store buffered write data in sample storage and reset the buffer
  for (auto &context : _accessor->writeDataContexts()) {
    context.storeBufferedData(time);
  }
}

void ParticipantImpl::trimOldDataBefore(double time)
{
  for (auto &context : _accessor->usedMeshContexts()) {
    for (const auto &name : context->mesh->availableData()) {
      context->mesh->data(name)->timeStepsStorage().trimBefore(time);
    }
  }
}

void ParticipantImpl::trimSendDataAfter(double time)
{
  for (auto &context : _accessor->writeDataContexts()) {
    context.trimAfter(time);
  }
}

void ParticipantImpl::finalize()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Finalized, "finalize() may only be called once.");

  // Events for the solver time, finally stopped here
  _solverAdvanceEvent.reset();

  Event e("finalize", profiling::Fundamental);

  if (_state == State::Initialized) {

    PRECICE_ASSERT(_couplingScheme->isInitialized());
    PRECICE_DEBUG("Finalize coupling scheme");
    _couplingScheme->finalize();

    closeCommunicationChannels(CloseChannels::All);
  }

  // Release ownership
  _couplingScheme.reset();
  _participants.clear();
  _accessor.reset();

  // Close Connections
  PRECICE_DEBUG("Close intra-participant communication");
  if (utils::IntraComm::isParallel()) {
    utils::IntraComm::getCommunication()->closeConnection();
    utils::IntraComm::getCommunication() = nullptr;
  }
  _m2ns.clear();

  // Stop and print Event logging
  e.stop();

  // Finalize PETSc and Events first
  utils::Petsc::finalize();
// This will lead to issues if we call finalize afterwards again
#ifndef PRECICE_NO_GINKGO
  device::Ginkgo::finalize();
#endif
  profiling::EventRegistry::instance().finalize();

  // Finally clear events and finalize MPI
  utils::Parallel::finalizeOrCleanupMPI();
  _state = State::Finalized;
}

int ParticipantImpl::getMeshDimensions(std::string_view meshName) const
{
  PRECICE_TRACE(meshName);
  PRECICE_VALIDATE_MESH_NAME(meshName);
  return _accessor->usedMeshContext(meshName).mesh->getDimensions();
}

int ParticipantImpl::getDataDimensions(std::string_view meshName, std::string_view dataName) const
{
  PRECICE_TRACE(meshName, dataName);
  PRECICE_VALIDATE_MESH_NAME(meshName);
  PRECICE_VALIDATE_DATA_NAME(meshName, dataName);
  return _accessor->usedMeshContext(meshName).mesh->data(dataName)->getDimensions();
}

bool ParticipantImpl::isCouplingOngoing() const
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Finalized, "isCouplingOngoing() cannot be called after finalize().");
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before isCouplingOngoing() can be evaluated.");
  return _couplingScheme->isCouplingOngoing();
}

bool ParticipantImpl::isTimeWindowComplete() const
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state != State::Constructed, "initialize() has to be called before isTimeWindowComplete().");
  PRECICE_CHECK(_state != State::Finalized, "isTimeWindowComplete() cannot be called after finalize().");
  return _couplingScheme->isTimeWindowComplete();
}

double ParticipantImpl::getMaxTimeStepSize() const
{
  PRECICE_CHECK(_state != State::Finalized, "getMaxTimeStepSize() cannot be called after finalize().");
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before getMaxTimeStepSize() can be evaluated.");
  const double nextTimeStepSize = _couplingScheme->getNextTimeStepMaxSize();
  // PRECICE_ASSERT(!math::equals(nextTimeStepSize, 0.0), nextTimeStepSize); // @todo requires https://github.com/precice/precice/issues/1904
  // PRECICE_ASSERT(math::greater(nextTimeStepSize, 0.0), nextTimeStepSize); // @todo requires https://github.com/precice/precice/issues/1904

  // safeguard needed because _couplingScheme->getNextTimeStepMaxSize() returns 0, if not isCouplingOngoing()
  // actual case where we want to warn the user
  PRECICE_WARN_IF(
      isCouplingOngoing() && not math::greater(nextTimeStepSize, 0.0, 100 * math::NUMERICAL_ZERO_DIFFERENCE),
      "preCICE just returned a maximum time step size of {}. Such a small value can happen if you use many substeps per time window over multiple time windows due to added-up differences of machine precision.",
      nextTimeStepSize);
  return nextTimeStepSize;
}

bool ParticipantImpl::requiresInitialData()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state == State::Constructed, "requiresInitialData() has to be called before initialize().");
  bool required = _couplingScheme->isActionRequired(cplscheme::CouplingScheme::Action::InitializeData);
  if (required) {
    _couplingScheme->markActionFulfilled(cplscheme::CouplingScheme::Action::InitializeData);
  }
  return required;
}

bool ParticipantImpl::requiresWritingCheckpoint()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before requiresWritingCheckpoint().");
  bool required = _couplingScheme->isActionRequired(cplscheme::CouplingScheme::Action::WriteCheckpoint);
  if (required) {
    _couplingScheme->markActionFulfilled(cplscheme::CouplingScheme::Action::WriteCheckpoint);
  }
  return required;
}

bool ParticipantImpl::requiresReadingCheckpoint()
{
  PRECICE_TRACE();
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before requiresReadingCheckpoint().");
  bool required = _couplingScheme->isActionRequired(cplscheme::CouplingScheme::Action::ReadCheckpoint);
  if (required) {
    _couplingScheme->markActionFulfilled(cplscheme::CouplingScheme::Action::ReadCheckpoint);
  }
  return required;
}

bool ParticipantImpl::requiresMeshConnectivityFor(std::string_view meshName) const
{
  PRECICE_VALIDATE_MESH_NAME(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  return context.meshRequirement == mapping::Mapping::MeshRequirement::FULL;
}

bool ParticipantImpl::requiresGradientDataFor(std::string_view meshName,
                                              std::string_view dataName) const
{
  PRECICE_VALIDATE_DATA_NAME(meshName, dataName);
  // Read data never requires gradients
  if (!_accessor->isDataWrite(meshName, dataName))
    return false;

  WriteDataContext &context = _accessor->writeDataContext(meshName, dataName);
  return context.hasGradient();
}

int ParticipantImpl::getMeshVertexSize(
    std::string_view meshName) const
{
  PRECICE_TRACE(meshName);
  PRECICE_REQUIRE_MESH_USE(meshName);
  // In case we access received mesh data: check, if the requested mesh data has already been received.
  // Otherwise, the function call doesn't make any sense
  PRECICE_CHECK((_state == State::Initialized) || _accessor->isMeshProvided(meshName), "initialize() has to be called before accessing"
                                                                                       " data of the received mesh \"{}\" on participant \"{}\".",
                meshName, _accessor->getName());
  MeshContext &context = _accessor->usedMeshContext(meshName);
  PRECICE_ASSERT(context.mesh.get() != nullptr);
  return context.mesh->nVertices();
}

/// @todo Currently not supported as we would need to re-compute the re-partition
void ParticipantImpl::resetMesh(
    std::string_view meshName)
{
  PRECICE_EXPERIMENTAL_API();
  PRECICE_CHECK(_allowsRemeshing, "Cannot reset meshes. This feature needs to be enabled using <precice-configuration experimental=\"1\" allow-remeshing=\"1\">.");
  PRECICE_CHECK(_state == State::Initialized, "initialize() has to be called before resetMesh().");
  PRECICE_TRACE(meshName);
  PRECICE_VALIDATE_MESH_NAME(meshName);
  PRECICE_CHECK(_couplingScheme->isCouplingOngoing(), "Cannot remesh after the last time window has been completed.");
  PRECICE_CHECK(_couplingScheme->isTimeWindowComplete(), "Cannot remesh while subcycling or iterating. Remeshing is only allowed when the time window is completed.");
  impl::MeshContext &context = _accessor->usedMeshContext(meshName);

  PRECICE_DEBUG("Clear mesh positions for mesh \"{}\"", context.mesh->getName());
  _meshLock.unlock(meshName);
  context.mesh->clear();
}

VertexID ParticipantImpl::setMeshVertex(
    std::string_view              meshName,
    ::precice::span<const double> position)
{
  PRECICE_TRACE(meshName);
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  auto &       mesh    = *context.mesh;
  PRECICE_CHECK(position.size() == static_cast<unsigned long>(mesh.getDimensions()),
                "Cannot set vertex for mesh \"{}\". Expected {} position components but found {}.", meshName, mesh.getDimensions(), position.size());
  Event e{fmt::format("setMeshVertex.{}", meshName), profiling::Fundamental};
  auto  index = mesh.createVertex(Eigen::Map<const Eigen::VectorXd>{position.data(), mesh.getDimensions()}).getID();
  mesh.allocateDataValues();

  const auto newSize = mesh.nVertices();
  for (auto &context : _accessor->writeDataContexts()) {
    if (context.getMeshName() == mesh.getName()) {
      context.resizeBufferTo(newSize);
    }
  }

  return index;
}

void ParticipantImpl::setMeshVertices(
    std::string_view              meshName,
    ::precice::span<const double> positions,
    ::precice::span<VertexID>     ids)
{
  PRECICE_TRACE(meshName, positions.size(), ids.size());
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  auto &       mesh    = *context.mesh;

  const auto meshDims             = mesh.getDimensions();
  const auto expectedPositionSize = ids.size() * meshDims;
  PRECICE_CHECK(positions.size() == expectedPositionSize,
                "Input sizes are inconsistent attempting to set vertices on {}D mesh \"{}\". "
                "You passed {} vertices indices and {} position components, but we expected {} position components ({} x {}).",
                meshDims, meshName, ids.size(), positions.size(), expectedPositionSize, ids.size(), meshDims);

  Event                                   e{fmt::format("setMeshVertices.{}", meshName), profiling::Fundamental};
  const Eigen::Map<const Eigen::MatrixXd> posMatrix{
      positions.data(), mesh.getDimensions(), static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(ids.size())};
  for (unsigned long i = 0; i < ids.size(); ++i) {
    ids[i] = mesh.createVertex(posMatrix.col(i)).getID();
  }
  mesh.allocateDataValues();

  const auto newSize = mesh.nVertices();
  for (auto &context : _accessor->writeDataContexts()) {
    if (context.getMeshName() == mesh.getName()) {
      context.resizeBufferTo(newSize);
    }
  }
}

void ParticipantImpl::setMeshEdge(
    std::string_view meshName,
    VertexID         first,
    VertexID         second)
{
  PRECICE_TRACE(meshName, first, second);
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh->isValidVertexID(first), errorInvalidVertexID(first));
    PRECICE_CHECK(mesh->isValidVertexID(second), errorInvalidVertexID(second));
    Event         e{fmt::format("setMeshEdge.{}", meshName), profiling::Fundamental};
    mesh::Vertex &v0 = mesh->vertex(first);
    mesh::Vertex &v1 = mesh->vertex(second);
    mesh->createEdge(v0, v1);
  }
}

void ParticipantImpl::setMeshEdges(
    std::string_view                meshName,
    ::precice::span<const VertexID> vertices)
{
  PRECICE_TRACE(meshName, vertices.size());
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement != mapping::Mapping::MeshRequirement::FULL) {
    return;
  }

  mesh::PtrMesh &mesh = context.mesh;
  PRECICE_CHECK(vertices.size() % 2 == 0,
                "Cannot interpret passed vertex IDs attempting to set edges of mesh \"{}\" . "
                "You passed {} vertices, but we expected an even number.",
                meshName, vertices.size());
  {
    auto end           = vertices.end();
    auto [first, last] = utils::find_first_range(vertices.begin(), end, [&mesh](VertexID vid) {
      return !mesh->isValidVertexID(vid);
    });
    PRECICE_CHECK(first == end,
                  impl::errorInvalidVertexIDRange,
                  std::distance(vertices.begin(), first),
                  std::distance(vertices.begin(), last));
  }

  Event e{fmt::format("setMeshEdges.{}", meshName), profiling::Fundamental};

  for (unsigned long i = 0; i < vertices.size() / 2; ++i) {
    auto aid = vertices[2 * i];
    auto bid = vertices[2 * i + 1];
    mesh->createEdge(mesh->vertex(aid), mesh->vertex(bid));
  }
}

void ParticipantImpl::setMeshTriangle(
    std::string_view meshName,
    VertexID         first,
    VertexID         second,
    VertexID         third)
{
  PRECICE_TRACE(meshName, first,
                second, third);

  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh->isValidVertexID(first), errorInvalidVertexID(first));
    PRECICE_CHECK(mesh->isValidVertexID(second), errorInvalidVertexID(second));
    PRECICE_CHECK(mesh->isValidVertexID(third), errorInvalidVertexID(third));
    PRECICE_CHECK(utils::unique_elements(utils::make_array(first, second, third)),
                  "setMeshTriangle() was called with repeated Vertex IDs ({}, {}, {}).",
                  first, second, third);
    mesh::Vertex *vertices[3];
    vertices[0] = &mesh->vertex(first);
    vertices[1] = &mesh->vertex(second);
    vertices[2] = &mesh->vertex(third);
    PRECICE_CHECK(utils::unique_elements(utils::make_array(vertices[0]->getCoords(),
                                                           vertices[1]->getCoords(), vertices[2]->getCoords())),
                  "setMeshTriangle() was called with vertices located at identical coordinates (IDs: {}, {}, {}).",
                  first, second, third);
    Event       e{fmt::format("setMeshTriangle.{}", meshName), profiling::Fundamental};
    mesh::Edge *edges[3];
    edges[0] = &mesh->createEdge(*vertices[0], *vertices[1]);
    edges[1] = &mesh->createEdge(*vertices[1], *vertices[2]);
    edges[2] = &mesh->createEdge(*vertices[2], *vertices[0]);

    mesh->createTriangle(*edges[0], *edges[1], *edges[2]);
  }
}

void ParticipantImpl::setMeshTriangles(
    std::string_view                meshName,
    ::precice::span<const VertexID> vertices)
{
  PRECICE_TRACE(meshName, vertices.size());
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement != mapping::Mapping::MeshRequirement::FULL) {
    return;
  }

  mesh::PtrMesh &mesh = context.mesh;
  PRECICE_CHECK(vertices.size() % 3 == 0,
                "Cannot interpret passed vertex IDs attempting to set triangles of mesh \"{}\" . "
                "You passed {} vertices, which isn't dividable by 3.",
                meshName, vertices.size());
  {
    auto end           = vertices.end();
    auto [first, last] = utils::find_first_range(vertices.begin(), end, [&mesh](VertexID vid) {
      return !mesh->isValidVertexID(vid);
    });
    PRECICE_CHECK(first == end,
                  impl::errorInvalidVertexIDRange,
                  std::distance(vertices.begin(), first),
                  std::distance(vertices.begin(), last));
  }

  Event e{fmt::format("setMeshTriangles.{}", meshName), profiling::Fundamental};

  for (unsigned long i = 0; i < vertices.size() / 3; ++i) {
    auto aid = vertices[3 * i];
    auto bid = vertices[3 * i + 1];
    auto cid = vertices[3 * i + 2];
    mesh->createTriangle(mesh->vertex(aid),
                         mesh->vertex(bid),
                         mesh->vertex(cid));
  }
}

void ParticipantImpl::setMeshQuad(
    std::string_view meshName,
    VertexID         first,
    VertexID         second,
    VertexID         third,
    VertexID         fourth)
{
  PRECICE_TRACE(meshName, first,
                second, third, fourth);
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  PRECICE_CHECK(context.mesh->getDimensions() == 3, "setMeshQuad is only possible for 3D meshes."
                                                    " Please set the mesh dimension to 3 in the preCICE configuration file.");
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    PRECICE_ASSERT(context.mesh);
    mesh::Mesh &mesh = *(context.mesh);
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh.isValidVertexID(first), errorInvalidVertexID(first));
    PRECICE_CHECK(mesh.isValidVertexID(second), errorInvalidVertexID(second));
    PRECICE_CHECK(mesh.isValidVertexID(third), errorInvalidVertexID(third));
    PRECICE_CHECK(mesh.isValidVertexID(fourth), errorInvalidVertexID(fourth));

    auto vertexIDs = utils::make_array(first, second, third, fourth);
    PRECICE_CHECK(utils::unique_elements(vertexIDs), "The four vertex ID's are not unique. Please check that the vertices that form the quad are correct.");

    auto coords = mesh::coordsFor(mesh, vertexIDs);
    PRECICE_CHECK(utils::unique_elements(coords),
                  "The four vertices that form the quad are not unique. The resulting shape may be a point, line or triangle."
                  "Please check that the adapter sends the four unique vertices that form the quad, or that the mesh on the interface is composed of quads.");

    auto convexity = math::geometry::isConvexQuad(coords);
    PRECICE_CHECK(convexity.convex, "The given quad is not convex. "
                                    "Please check that the adapter send the four correct vertices or that the interface is composed of quads.");
    auto reordered = utils::reorder_array(convexity.vertexOrder, mesh::vertexPtrsFor(mesh, vertexIDs));

    Event e{fmt::format("setMeshQuad.{}", meshName), profiling::Fundamental};

    // Vertices are now in the order: V0-V1-V2-V3-V0.
    // Use the shortest diagonal to split the quad into 2 triangles.
    // Vertices are now in V0-V1-V2-V3-V0 order. The new edge, e[4] is either 0-2 or 1-3
    double distance02 = (reordered[0]->getCoords() - reordered[2]->getCoords()).norm();
    double distance13 = (reordered[1]->getCoords() - reordered[3]->getCoords()).norm();

    // The new edge, e[4], is the shortest diagonal of the quad
    if (distance02 <= distance13) {
      mesh.createTriangle(*reordered[0], *reordered[2], *reordered[1]);
      mesh.createTriangle(*reordered[0], *reordered[2], *reordered[3]);
    } else {
      mesh.createTriangle(*reordered[1], *reordered[3], *reordered[0]);
      mesh.createTriangle(*reordered[1], *reordered[3], *reordered[2]);
    }
  }
}

void ParticipantImpl::setMeshQuads(
    std::string_view                meshName,
    ::precice::span<const VertexID> vertices)
{
  PRECICE_TRACE(meshName, vertices.size());
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement != mapping::Mapping::MeshRequirement::FULL) {
    return;
  }

  mesh::Mesh &mesh = *(context.mesh);
  PRECICE_CHECK(vertices.size() % 4 == 0,
                "Cannot interpret passed vertex IDs attempting to set quads of mesh \"{}\" . "
                "You passed {} vertices, which isn't dividable by 4.",
                meshName, vertices.size());
  {
    auto end           = vertices.end();
    auto [first, last] = utils::find_first_range(vertices.begin(), end, [&mesh](VertexID vid) {
      return !mesh.isValidVertexID(vid);
    });
    PRECICE_CHECK(first == end,
                  impl::errorInvalidVertexIDRange,
                  std::distance(vertices.begin(), first),
                  std::distance(vertices.begin(), last));
  }

  for (unsigned long i = 0; i < vertices.size() / 4; ++i) {
    auto aid = vertices[4 * i];
    auto bid = vertices[4 * i + 1];
    auto cid = vertices[4 * i + 2];
    auto did = vertices[4 * i + 3];

    auto vertexIDs = utils::make_array(aid, bid, cid, did);
    PRECICE_CHECK(utils::unique_elements(vertexIDs), "The four vertex ID's of the quad nr {} are not unique. Please check that the vertices that form the quad are correct.", i);

    auto coords = mesh::coordsFor(mesh, vertexIDs);
    PRECICE_CHECK(utils::unique_elements(coords),
                  "The four vertices that form the quad nr {} are not unique. The resulting shape may be a point, line or triangle."
                  "Please check that the adapter sends the four unique vertices that form the quad, or that the mesh on the interface is composed of quads.",
                  i);

    auto convexity = math::geometry::isConvexQuad(coords);
    PRECICE_CHECK(convexity.convex, "The given quad nr {} is not convex. "
                                    "Please check that the adapter send the four correct vertices or that the interface is composed of quads.",
                  i);
    auto reordered = utils::reorder_array(convexity.vertexOrder, mesh::vertexPtrsFor(mesh, vertexIDs));

    Event e{fmt::format("setMeshQuads.{}", meshName), profiling::Fundamental};

    // Use the shortest diagonal to split the quad into 2 triangles.
    // Vertices are now in V0-V1-V2-V3-V0 order. The new edge, e[4] is either 0-2 or 1-3
    double distance02 = (reordered[0]->getCoords() - reordered[2]->getCoords()).norm();
    double distance13 = (reordered[1]->getCoords() - reordered[3]->getCoords()).norm();

    if (distance02 <= distance13) {
      mesh.createTriangle(*reordered[0], *reordered[2], *reordered[1]);
      mesh.createTriangle(*reordered[0], *reordered[2], *reordered[3]);
    } else {
      mesh.createTriangle(*reordered[1], *reordered[3], *reordered[0]);
      mesh.createTriangle(*reordered[1], *reordered[3], *reordered[2]);
    }
  }
}

void ParticipantImpl::setMeshTetrahedron(
    std::string_view meshName,
    VertexID         first,
    VertexID         second,
    VertexID         third,
    VertexID         fourth)
{
  PRECICE_TRACE(meshName, first, second, third, fourth);
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  PRECICE_CHECK(context.mesh->getDimensions() == 3, "setMeshTetrahedron is only possible for 3D meshes."
                                                    " Please set the mesh dimension to 3 in the preCICE configuration file.");
  Event e{fmt::format("setMeshTetrahedron.{}", meshName), profiling::Fundamental};

  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL) {
    mesh::PtrMesh &mesh = context.mesh;
    using impl::errorInvalidVertexID;
    PRECICE_CHECK(mesh->isValidVertexID(first), errorInvalidVertexID(first));
    PRECICE_CHECK(mesh->isValidVertexID(second), errorInvalidVertexID(second));
    PRECICE_CHECK(mesh->isValidVertexID(third), errorInvalidVertexID(third));
    PRECICE_CHECK(mesh->isValidVertexID(fourth), errorInvalidVertexID(fourth));
    mesh::Vertex &A = mesh->vertex(first);
    mesh::Vertex &B = mesh->vertex(second);
    mesh::Vertex &C = mesh->vertex(third);
    mesh::Vertex &D = mesh->vertex(fourth);

    mesh->createTetrahedron(A, B, C, D);
  }
}

void ParticipantImpl::setMeshTetrahedra(
    std::string_view                meshName,
    ::precice::span<const VertexID> vertices)
{
  PRECICE_TRACE(meshName, vertices.size());
  PRECICE_REQUIRE_MESH_MODIFY(meshName);
  MeshContext &context = _accessor->usedMeshContext(meshName);
  if (context.meshRequirement != mapping::Mapping::MeshRequirement::FULL) {
    return;
  }

  mesh::PtrMesh &mesh = context.mesh;
  PRECICE_CHECK(vertices.size() % 4 == 0,
                "Cannot interpret passed vertex IDs attempting to set quads of mesh \"{}\" . "
                "You passed {} vertices, which isn't dividable by 4.",
                meshName, vertices.size());
  {
    auto end           = vertices.end();
    auto [first, last] = utils::find_first_range(vertices.begin(), end, [&mesh](VertexID vid) {
      return !mesh->isValidVertexID(vid);
    });
    PRECICE_CHECK(first == end,
                  impl::errorInvalidVertexIDRange,
                  std::distance(vertices.begin(), first),
                  std::distance(vertices.begin(), last));
  }

  Event e{fmt::format("setMeshTetrahedra.{}", meshName), profiling::Fundamental};

  for (unsigned long i = 0; i < vertices.size() / 4; ++i) {
    auto aid = vertices[4 * i];
    auto bid = vertices[4 * i + 1];
    auto cid = vertices[4 * i + 2];
    auto did = vertices[4 * i + 3];
    mesh->createTetrahedron(mesh->vertex(aid),
                            mesh->vertex(bid),
                            mesh->vertex(cid),
                            mesh->vertex(did));
  }
}

void ParticipantImpl::writeData(
    std::string_view                meshName,
    std::string_view                dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   values)
{
  PRECICE_TRACE(meshName, dataName, vertices.size());
  PRECICE_CHECK(_state != State::Finalized, "writeData(...) cannot be called after finalize().");
  PRECICE_CHECK(_state == State::Constructed || (_state == State::Initialized && isCouplingOngoing()), "Calling writeData(...) is forbidden if coupling is not ongoing, because the data you are trying to write will not be used anymore. You can fix this by always calling writeData(...) before the advance(...) call in your simulation loop or by using Participant::isCouplingOngoing() to implement a safeguard.");
  PRECICE_REQUIRE_DATA_WRITE(meshName, dataName);
  // Inconsistent sizes will be handled below
  if (vertices.empty() && values.empty()) {
    return;
  }

  WriteDataContext &context = _accessor->writeDataContext(meshName, dataName);

  const auto dataDims         = context.getDataDimensions();
  const auto expectedDataSize = vertices.size() * dataDims;
  PRECICE_CHECK(expectedDataSize == values.size(),
                "Input sizes are inconsistent attempting to write {}D data \"{}\" to mesh \"{}\". "
                "You passed {} vertices and {} data components, but we expected {} data components ({} x {}).",
                dataDims, dataName, meshName,
                vertices.size(), values.size(), expectedDataSize, dataDims, vertices.size());

  // Sizes are correct at this point
  PRECICE_VALIDATE_DATA(values.data(), values.size()); // TODO Only take span

  if (auto index = context.locateInvalidVertexID(vertices); index) {
    PRECICE_ERROR("Cannot write data \"{}\" to mesh \"{}\" due to invalid Vertex ID at vertices[{}]. "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  dataName, meshName, *index);
  }
  Event e{fmt::format("writeData.{}_{}", meshName, dataName), profiling::Fundamental};
  context.writeValuesIntoDataBuffer(vertices, values);
}

void ParticipantImpl::readData(
    std::string_view                meshName,
    std::string_view                dataName,
    ::precice::span<const VertexID> vertices,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  PRECICE_TRACE(meshName, dataName, vertices.size(), relativeReadTime);
  PRECICE_CHECK(_state != State::Constructed, "readData(...) cannot be called before initialize().");
  PRECICE_CHECK(_state != State::Finalized, "readData(...) cannot be called after finalize().");
  PRECICE_CHECK(math::smallerEquals(relativeReadTime, _couplingScheme->getNextTimeStepMaxSize()), "readData(...) cannot sample data outside of current time window.");
  PRECICE_CHECK(relativeReadTime >= 0, "readData(...) cannot sample data before the current time.");
  PRECICE_CHECK(isCouplingOngoing() || math::equals(relativeReadTime, 0.0), "Calling readData(...) with relativeReadTime = {} is forbidden if coupling is not ongoing. If coupling finished, only data for relativeReadTime = 0 is available. Please always use precice.getMaxTimeStepSize() to obtain the maximum allowed relativeReadTime.", relativeReadTime);

  PRECICE_REQUIRE_DATA_READ(meshName, dataName);

  PRECICE_CHECK(_meshLock.check(meshName),
                "Cannot read from mesh \"{}\" after it has been reset. Please read data before calling resetMesh().",
                meshName);

  // Inconsistent sizes will be handled below
  if (vertices.empty() && values.empty()) {
    return;
  }

  ReadDataContext &context = _accessor->readDataContext(meshName, dataName);
  PRECICE_CHECK(context.hasSamples(), "Data \"{}\" cannot be read from mesh \"{}\" as it contains no samples. "
                                      "This is typically a configuration issue of the data flow. "
                                      "Check if the data is correctly exchanged to this participant \"{}\" and mapped to mesh \"{}\".",
                dataName, meshName, _accessorName, meshName);

  const auto dataDims         = context.getDataDimensions();
  const auto expectedDataSize = vertices.size() * dataDims;
  PRECICE_CHECK(expectedDataSize == values.size(),
                "Input/Output sizes are inconsistent attempting to read {}D data \"{}\" from mesh \"{}\". "
                "You passed {} vertices and {} data components, but we expected {} data components ({} x {}).",
                dataDims, dataName, meshName,
                vertices.size(), values.size(), expectedDataSize, dataDims, vertices.size());

  if (auto index = context.locateInvalidVertexID(vertices); index) {
    PRECICE_ERROR("Cannot read data \"{}\" from mesh \"{}\" due to invalid Vertex ID at vertices[{}]. "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  dataName, meshName, *index);
  }

  Event e{fmt::format("readData.{}_{}", meshName, dataName), profiling::Fundamental};

  double readTime = _couplingScheme->getTime() + relativeReadTime;
  context.readValues(vertices, readTime, values);
}

//////////////////////////////////////////////
void ParticipantImpl::mapAndreadData(
    std::string_view              meshName,
    std::string_view              dataName,
    ::precice::span<const double> coordinates,
    double                        relativeReadTime,
    ::precice::span<double>       values) const
{
  PRECICE_TRACE(meshName, dataName, coordinates.size(), relativeReadTime);
  // TODO: Make these checks conditional
  PRECICE_CHECK(_state != State::Constructed, "mapAndreadData(...) cannot be called before initialize().");
  PRECICE_CHECK(_state != State::Finalized, "mapAndreadData(...) cannot be called after finalize().");
  PRECICE_CHECK(math::smallerEquals(relativeReadTime, _couplingScheme->getNextTimeStepMaxSize()), "readData(...) cannot sample data outside of current time window.");
  PRECICE_CHECK(relativeReadTime >= 0, "mapAndreadData(...) cannot sample data before the current time.");
  PRECICE_CHECK(isCouplingOngoing() || math::equals(relativeReadTime, 0.0), "Calling mapAndreadData(...) with relativeReadTime = {} is forbidden if coupling is not ongoing. If coupling finished, only data for relativeReadTime = 0 is available. Please always use precice.getMaxTimeStepSize() to obtain the maximum allowed relativeReadTime.", relativeReadTime);

  PRECICE_REQUIRE_DATA_READ(meshName, dataName);
  PRECICE_VALIDATE_DATA(coordinates.begin(), coordinates.size());

  // Inconsistent sizes will be handled below
  if (coordinates.empty() && values.empty()) {
    return;
  }

  // Make use of the read data context
  // Note that meshName refers here typically to a remote mesh
  ReadDataContext &dataContext      = _accessor->readDataContext(meshName, dataName);
  const auto       dataDims         = dataContext.getDataDimensions();
  const auto       dim              = dataContext.getSpatialDimensions();
  const auto       expectedDataSize = (coordinates.size() / dim) * dataDims;
  // TODO: Add check that this vertex is within the access region?
  double readTime = _couplingScheme->getTime() + relativeReadTime;
  dataContext.mapAndReadValues(coordinates, readTime, values);
}
////////////////////////////////////////////////////

void ParticipantImpl::mapAndwriteData(
    std::string_view              meshName,
    std::string_view              dataName,
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  PRECICE_TRACE(meshName, dataName, coordinates.size());
  PRECICE_CHECK(_state != State::Finalized, "mapAndwriteData(...) cannot be called after finalize().");
  PRECICE_CHECK(_state == State::Constructed || (_state == State::Initialized && isCouplingOngoing()), "Calling mapAndwriteData(...) is forbidden if coupling is not ongoing, because the data you are trying to write will not be used anymore. You can fix this by always calling mapAndwriteData(...) before the advance(...) call in your simulation loop or by using Participant::isCouplingOngoing() to implement a safeguard.");
  PRECICE_REQUIRE_DATA_WRITE(meshName, dataName);
  // Inconsistent sizes will be handled below
  if (coordinates.empty() && values.empty()) {
    return;
  }

  WriteDataContext &dataContext = _accessor->writeDataContext(meshName, dataName);

  const auto dataDims         = dataContext.getDataDimensions();
  const auto dim              = dataContext.getSpatialDimensions();
  const auto expectedDataSize = (coordinates.size() / dim) * dataDims;
  PRECICE_CHECK(expectedDataSize == values.size(),
                "Input sizes are inconsistent attempting to write {}D data \"{}\" to mesh \"{}\". "
                "You passed {} vertices and {} data components, but we expected {} data components ({} x {}).",
                dataDims, dataName, meshName,
                coordinates.size() / dim, values.size(), expectedDataSize, dataDims, coordinates.size() / dim);

  // Sizes are correct at this point
  PRECICE_VALIDATE_DATA(values.data(), values.size()); // TODO Only take span

  // TODO: Add check that this vertex is within the access region?
  dataContext.mapAndWriteValues(coordinates, values);
}

void ParticipantImpl::writeGradientData(
    std::string_view                meshName,
    std::string_view                dataName,
    ::precice::span<const VertexID> vertices,
    ::precice::span<const double>   gradients)
{
  PRECICE_EXPERIMENTAL_API();

  // Asserts and checks
  PRECICE_TRACE(meshName, dataName, vertices.size());
  PRECICE_CHECK(_state != State::Finalized, "writeGradientData(...) cannot be called after finalize().");
  PRECICE_REQUIRE_DATA_WRITE(meshName, dataName);

  // Inconsistent sizes will be handled below
  if ((vertices.empty() && gradients.empty()) || !requiresGradientDataFor(meshName, dataName)) {
    return;
  }

  // Get the data
  WriteDataContext &context = _accessor->writeDataContext(meshName, dataName);

  // Check if the Data object of given mesh has been initialized with gradient data
  PRECICE_CHECK(context.hasGradient(), "Data \"{}\" has no gradient values available. Please set the gradient flag to true under the data attribute in the configuration file.", dataName);

  if (auto index = context.locateInvalidVertexID(vertices); index) {
    PRECICE_ERROR("Cannot write gradient data \"{}\" to mesh \"{}\" due to invalid Vertex ID at vertices[{}]. "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  dataName, meshName, *index);
  }

  const auto dataDims           = context.getDataDimensions();
  const auto meshDims           = context.getSpatialDimensions();
  const auto gradientComponents = meshDims * dataDims;
  const auto expectedComponents = vertices.size() * gradientComponents;
  PRECICE_CHECK(expectedComponents == gradients.size(),
                "Input sizes are inconsistent attempting to write gradient for data \"{}\" to mesh \"{}\". "
                "A single gradient/Jacobian for {}D data on a {}D mesh has {} components. "
                "You passed {} vertices and {} gradient components, but we expected {} gradient components. ",
                dataName, meshName,
                dataDims, meshDims, gradientComponents,
                vertices.size(), gradients.size(), expectedComponents);

  PRECICE_VALIDATE_DATA(gradients.data(), gradients.size());

  Event e{fmt::format("writeGradientData.{}_{}", meshName, dataName), profiling::Fundamental};

  context.writeGradientsIntoDataBuffer(vertices, gradients);
}

void ParticipantImpl::setMeshAccessRegion(
    const std::string_view        meshName,
    ::precice::span<const double> boundingBox) const
{
  PRECICE_TRACE(meshName, boundingBox.size());
  PRECICE_REQUIRE_MESH_USE(meshName);
  PRECICE_CHECK(_state != State::Finalized, "setMeshAccessRegion() cannot be called after finalize().");
  PRECICE_CHECK(_state != State::Initialized, "setMeshAccessRegion() needs to be called before initialize().");

  // Get the related mesh
  MeshContext &context = _accessor->meshContext(meshName);

  PRECICE_CHECK(!context.accessRegionDefined, "A mesh access region was already defined for mesh \"{}\". setMeshAccessRegion may only be called once per mesh.", context.mesh->getName());
  mesh::PtrMesh mesh(context.mesh);
  int           dim = mesh->getDimensions();
  PRECICE_CHECK(boundingBox.size() == static_cast<unsigned long>(dim) * 2,
                "Incorrect amount of bounding box components attempting to set the bounding box of {}D mesh \"{}\" . "
                "You passed {} limits, but we expected {} ({}x2).",
                dim, meshName, boundingBox.size(), dim * 2, dim);

  // Transform bounds into a suitable format
  PRECICE_DEBUG("Define bounding box");
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
  context.accessRegionDefined = true;
}

void ParticipantImpl::getMeshVertexIDsAndCoordinates(
    const std::string_view    meshName,
    ::precice::span<VertexID> ids,
    ::precice::span<double>   coordinates) const
{
  PRECICE_TRACE(meshName, ids.size(), coordinates.size());
  PRECICE_REQUIRE_MESH_USE(meshName);
  PRECICE_DEBUG("Get {} mesh vertices with IDs", ids.size());

  // Check, if the requested mesh data has already been received. Otherwise, the function call doesn't make any sense
  PRECICE_CHECK((_state == State::Initialized) || _accessor->isMeshProvided(meshName), "initialize() has to be called before accessing"
                                                                                       " data of the received mesh \"{}\" on participant \"{}\".",
                meshName, _accessor->getName());

  if (ids.empty() && coordinates.empty()) {
    return;
  }
  const MeshContext & context = _accessor->meshContext(meshName);
  const mesh::PtrMesh mesh(context.mesh);

  const auto meshSize = mesh->nVertices();
  const auto meshDims = mesh->getDimensions();
  PRECICE_CHECK(ids.size() == meshSize,
                "Output size is incorrect attempting to get vertex ids of {}D mesh \"{}\". "
                "You passed {} vertices indices, but we expected {}. "
                "Use getMeshVertexSize(\"{}\") to receive the required amount of vertices.",
                meshDims, meshName, ids.size(), meshSize, meshName);
  const auto expectedCoordinatesSize = static_cast<unsigned long>(meshDims * meshSize);
  PRECICE_CHECK(coordinates.size() == expectedCoordinatesSize,
                "Output size is incorrect attempting to get vertex coordinates of {}D mesh \"{}\". "
                "You passed {} coordinate components, but we expected {} ({}x{}). "
                "Use getMeshVertexSize(\"{}\") and getMeshDimensions(\"{}\") to receive the required amount components",
                meshDims, meshName, coordinates.size(), expectedCoordinatesSize, meshSize, meshDims, meshName, meshName);

  PRECICE_CHECK(ids.size() <= meshSize, "The queried size exceeds the number of available points.");
  Event e{fmt::format("getMeshVertexIDsAndCoordinates.{}", meshName), profiling::Fundamental};

  Eigen::Map<Eigen::MatrixXd> posMatrix{
      coordinates.data(), mesh->getDimensions(), static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(ids.size())};

  for (unsigned long i = 0; i < ids.size(); i++) {
    PRECICE_ASSERT(mesh->isValidVertexID(i), i, meshSize);
    ids[i]           = mesh->vertex(i).getID();
    posMatrix.col(i) = mesh->vertex(i).getCoords();
  }
}

void ParticipantImpl::configureM2Ns(
    const m2n::M2NConfiguration::SharedPointer &config)
{
  PRECICE_TRACE();
  for (const auto &m2nConf : config->m2ns()) {
    if (m2nConf.acceptor != _accessorName && m2nConf.connector != _accessorName) {
      continue;
    }

    std::string comPartner("");
    bool        isRequesting;
    if (m2nConf.acceptor == _accessorName) {
      comPartner   = m2nConf.connector;
      isRequesting = true;
    } else {
      comPartner   = m2nConf.acceptor;
      isRequesting = false;
    }

    PRECICE_ASSERT(!comPartner.empty());
    for (const impl::PtrParticipant &participant : _participants) {
      if (participant->getName() == comPartner) {
        PRECICE_ASSERT(not utils::contained(comPartner, _m2ns), comPartner);
        PRECICE_ASSERT(m2nConf.m2n);

        _m2ns[comPartner] = [&] {
          m2n::BoundM2N bound;
          bound.m2n          = m2nConf.m2n;
          bound.localName    = _accessorName;
          bound.remoteName   = comPartner;
          bound.isRequesting = isRequesting;
          return bound;
        }();
      }
    }
  }
}

void ParticipantImpl::configurePartitions(
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

    } else { // Accessor receives mesh
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

void ParticipantImpl::compareBoundingBoxes()
{
  // sort meshContexts by name, for communication in right order.
  std::sort(_accessor->usedMeshContexts().begin(), _accessor->usedMeshContexts().end(),
            [](MeshContext const *const lhs, MeshContext const *const rhs) -> bool {
              return lhs->mesh->getName() < rhs->mesh->getName();
            });

  for (MeshContext *meshContext : _accessor->usedMeshContexts()) {
    if (meshContext->provideMesh) // provided meshes need their bounding boxes already for the re-partitioning
      meshContext->mesh->computeBoundingBox();

    meshContext->clearMappings();
  }

  for (MeshContext *meshContext : _accessor->usedMeshContexts()) {
    meshContext->partition->compareBoundingBoxes();
  }
}

void ParticipantImpl::computePartitions()
{
  // We need to do this in two loops: First, communicate the mesh and later compute the partition.
  // Originally, this was done in one loop. This however gave deadlock if two meshes needed to be communicated cross-wise.
  // Both loops need a different sorting

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

    const auto requiredSize = meshContext->mesh->nVertices();
    for (auto &context : _accessor->writeDataContexts()) {
      if (context.getMeshName() == meshContext->mesh->getName()) {
        context.resizeBufferTo(requiredSize);
      }
    }
  }
}

void ParticipantImpl::computeMappings(std::vector<MappingContext> &contexts, const std::string &mappingType)
{
  PRECICE_TRACE();
  using namespace mapping;
  for (impl::MappingContext &context : contexts) {
    if (not context.mapping->hasComputedMapping()) {
      PRECICE_INFO_IF(context.configuredWithAliasTag,
                      "Automatic RBF mapping alias from mesh \"{}\" to mesh \"{}\" in \"{}\" direction resolves to \"{}\" .",
                      context.mapping->getInputMesh()->getName(), context.mapping->getOutputMesh()->getName(), mappingType, context.mapping->getName());
      PRECICE_INFO("Computing \"{}\" mapping from mesh \"{}\" to mesh \"{}\" in \"{}\" direction.",
                   context.mapping->getName(), context.mapping->getInputMesh()->getName(), context.mapping->getOutputMesh()->getName(), mappingType);
      context.mapping->computeMapping();
    }
  }
}

void ParticipantImpl::mapInitialWrittenData()
{
  PRECICE_TRACE();
  if (!_accessor->hasWriteMappings()) {
    return;
  }

  Event e("mapping", profiling::Fundamental);
  computeMappings(_accessor->writeMappingContexts(), "write");
  for (auto &context : _accessor->writeDataContexts()) {
    if (context.hasMapping()) {
      PRECICE_DEBUG("Map initial write data \"{}\" from mesh \"{}\"", context.getDataName(), context.getMeshName());
      _executedWriteMappings += context.mapData(std::nullopt, true);
    }
  }
}

void ParticipantImpl::mapWrittenData(std::optional<double> after)
{
  PRECICE_TRACE();
  if (!_accessor->hasWriteMappings()) {
    return;
  }

  Event e("mapping", profiling::Fundamental);
  computeMappings(_accessor->writeMappingContexts(), "write");
  for (auto &context : _accessor->writeDataContexts()) {
    if (context.hasMapping()) {
      PRECICE_DEBUG("Map write data \"{}\" from mesh \"{}\"", context.getDataName(), context.getMeshName());
      _executedWriteMappings += context.mapData(after);
    }
  }
}

void ParticipantImpl::trimReadMappedData(double startOfTimeWindow, bool isTimeWindowComplete, const cplscheme::ImplicitData &fromData)
{
  PRECICE_TRACE();
  for (auto &context : _accessor->readDataContexts()) {
    if (context.hasMapping()) {
      if (isTimeWindowComplete) {
        // For serial implicit second, we need to discard everything before startOfTimeWindow to preserve the time window start
        // For serial implicit first, we need to discard everything as everything is new
        // For parallel implicit, we need to discard everything as everything is new
        context.clearToDataFor(fromData);
      } else {
        context.trimToDataAfterFor(fromData, startOfTimeWindow);
      }
    }
  }
}

void ParticipantImpl::mapInitialReadData()
{
  PRECICE_TRACE();
  if (!_accessor->hasReadMappings()) {
    return;
  }

  Event e("mapping", profiling::Fundamental);
  computeMappings(_accessor->readMappingContexts(), "read");
  for (auto &context : _accessor->readDataContexts()) {
    if (context.hasMapping()) {
      PRECICE_DEBUG("Map initial read data \"{}\" to mesh \"{}\"", context.getDataName(), context.getMeshName());
      // We always ensure that all read data was mapped
      _executedReadMappings += context.mapData(std::nullopt, true);
    }
  }
}

void ParticipantImpl::mapReadData()
{
  PRECICE_TRACE();
  if (!_accessor->hasReadMappings()) {
    return;
  }

  Event e("mapping", profiling::Fundamental);
  computeMappings(_accessor->readMappingContexts(), "read");
  for (auto &context : _accessor->readDataContexts()) {
    if (context.hasMapping()) {
      PRECICE_DEBUG("Map read data \"{}\" to mesh \"{}\"", context.getDataName(), context.getMeshName());
      // We always ensure that all read data was mapped
      _executedReadMappings += context.mapData();
    }
  }
}

void ParticipantImpl::performDataActions(const std::set<action::Action::Timing> &timings)
{
  PRECICE_TRACE();
  for (action::PtrAction &action : _accessor->actions()) {
    if (timings.find(action->getTiming()) != timings.end()) {
      action->performAction();
    }
  }
}

void ParticipantImpl::handleExports(ExportTiming timing)
{
  PRECICE_TRACE();
  if (!_accessor->hasExports()) {
    return;
  }
  PRECICE_DEBUG("Handle exports");
  profiling::Event e{"handleExports"};

  if (timing == ExportTiming::Initial) {
    _accessor->exportInitial();
    return;
  }

  ParticipantState::IntermediateExport exp;
  exp.timewindow = _couplingScheme->getTimeWindows() - 1;
  exp.iteration  = _numberAdvanceCalls;
  exp.complete   = _couplingScheme->isTimeWindowComplete();
  exp.final      = !_couplingScheme->isCouplingOngoing();
  exp.time       = _couplingScheme->getTime();
  _accessor->exportIntermediate(exp);
}

void ParticipantImpl::resetWrittenData()
{
  PRECICE_TRACE();
  for (auto &context : _accessor->writeDataContexts()) {
    context.resetBuffer();
  }
}

PtrParticipant ParticipantImpl::determineAccessingParticipant(
    const config::Configuration &config)
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

void ParticipantImpl::initializeIntraCommunication()
{
  PRECICE_TRACE();

  Event e("com.initializeIntraCom", profiling::Fundamental);
  utils::IntraComm::getCommunication()->connectIntraComm(
      _accessorName, "IntraComm",
      _accessorProcessRank, _accessorCommunicatorSize);
  utils::IntraComm::barrier();
}

void ParticipantImpl::syncTimestep(double computedTimeStepSize)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel());
  Event e("syncTimestep", profiling::Fundamental);
  if (utils::IntraComm::isSecondary()) {
    utils::IntraComm::getCommunication()->send(computedTimeStepSize, 0);
  } else {
    PRECICE_ASSERT(utils::IntraComm::isPrimary());
    for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      double dt;
      utils::IntraComm::getCommunication()->receive(dt, secondaryRank);
      PRECICE_CHECK(math::equals(dt, computedTimeStepSize),
                    "Found ambiguous values for the time step size passed to preCICE in \"advance\". On rank {}, the value is {}, while on rank 0, the value is {}.",
                    secondaryRank, dt, computedTimeStepSize);
    }
  }
}

void ParticipantImpl::advanceCouplingScheme()
{
  PRECICE_DEBUG("Advance coupling scheme");
  Event e("advanceCoupling", profiling::Fundamental);
  // Orchestrate local and remote mesh changes
  std::vector<MeshID> localChanges;

  [[maybe_unused]] auto remoteChanges1 = _couplingScheme->firstSynchronization(localChanges);
  _couplingScheme->firstExchange();
  // Orchestrate remote mesh changes (local ones were handled in the first sync)
  [[maybe_unused]] auto remoteChanges2 = _couplingScheme->secondSynchronization();
  _couplingScheme->secondExchange();
}

void ParticipantImpl::closeCommunicationChannels(CloseChannels close)
{
  // Optionally apply some final ping-pong to sync solver that run e.g. with a uni-directional coupling
  // afterwards close connections
  PRECICE_INFO("{} {}communication channels",
               (_waitInFinalize ? "Synchronize participants and close" : "Close"),
               (close == CloseChannels::Distributed ? "distributed " : ""));
  std::string ping = "ping";
  std::string pong = "pong";
  for (auto &iter : _m2ns) {
    auto bm2n = iter.second;
    if (!bm2n.m2n->isConnected()) {
      PRECICE_DEBUG("Skipping closure of defective connection with {}", bm2n.remoteName);
      continue;
    }
    if (_waitInFinalize && not utils::IntraComm::isSecondary()) {
      auto comm = bm2n.m2n->getPrimaryRankCommunication();
      PRECICE_DEBUG("Synchronizing primary rank with {}", bm2n.remoteName);
      if (bm2n.isRequesting) {
        comm->send(ping, 0);
        std::string receive = "init";
        comm->receive(receive, 0);
        PRECICE_ASSERT(receive == pong);
      } else {
        std::string receive = "init";
        comm->receive(receive, 0);
        PRECICE_ASSERT(receive == ping);
        comm->send(pong, 0);
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

const mesh::Mesh &ParticipantImpl::mesh(const std::string &meshName) const
{
  PRECICE_TRACE(meshName);
  return *_accessor->usedMeshContext(meshName).mesh;
}

ParticipantImpl::MappedSamples ParticipantImpl::mappedSamples() const
{
  MappedSamples res;
  res.read  = _executedReadMappings;
  res.write = _executedWriteMappings;
  return res;
}

// Reinitialization

int ParticipantImpl::getTotalMeshChanges() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_allowsRemeshing);
  Event e("remesh.exchangeLocalMeshChanges", profiling::Synchronize);
  int   localMeshesChanges = _meshLock.countUnlocked();
  PRECICE_DEBUG("Mesh changes of rank: {}", localMeshesChanges);

  int totalMeshesChanges = 0;
  utils::IntraComm::allreduceSum(localMeshesChanges, totalMeshesChanges);
  PRECICE_DEBUG("Mesh changes of participant: {}", totalMeshesChanges);
  return totalMeshesChanges;
}

bool ParticipantImpl::reinitHandshake(bool requestReinit) const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_allowsRemeshing);
  Event e("remesh.exchangeRemoteMeshChanges", profiling::Synchronize);

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Remeshing is{} required by this participant.", (requestReinit ? "" : " not"));

    bool swarmReinitRequired = requestReinit;
    for (auto &iter : _m2ns) {
      PRECICE_DEBUG("Coordinating remeshing with {}", iter.first);
      bool  received = false;
      auto &comm     = *iter.second.m2n->getPrimaryRankCommunication();
      if (iter.second.isRequesting) {
        comm.send(requestReinit, 0);
        comm.receive(received, 0);
      } else {
        comm.receive(received, 0);
        comm.send(requestReinit, 0);
      }
      swarmReinitRequired |= received;
    }
    PRECICE_DEBUG("Coordinated that overall{} remeshing is required.", (swarmReinitRequired ? "" : " no"));

    utils::IntraComm::broadcast(swarmReinitRequired);
    return swarmReinitRequired;
  } else {
    bool swarmReinitRequired = false;
    utils::IntraComm::broadcast(swarmReinitRequired);
    return swarmReinitRequired;
  }
}

} // namespace precice::impl

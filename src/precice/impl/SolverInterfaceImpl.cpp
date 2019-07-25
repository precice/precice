#include "SolverInterfaceImpl.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "precice/impl/RequestManager.hpp"
#include "precice/impl/ValidationMacros.hpp"
#include "precice/config/Configuration.hpp"
#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Merge.hpp"
#include "io/ExportContext.hpp"
#include "io/Export.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/BoundM2N.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "utils/EventUtils.hpp"
#include "utils/Helpers.hpp"
#include "utils/SignalHandler.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "mapping/Mapping.hpp"
#include <Eigen/Core>
#include "partition/ReceivedPartition.hpp"
#include "partition/ProvidedPartition.hpp"
#include "versions.hpp"

#include <csignal> // used for installing crash handler
#ifndef SIGXCPU
#define SIGXCPU 24 /* exceeded CPU time limit */
#endif

#include <utility>
#include <algorithm>

#include "logging/Logger.hpp"
#include "logging/LogConfiguration.hpp"


using precice::utils::Event;
using precice::utils::EventRegistry;

namespace precice {

/// Set to true if unit/integration tests are executed
bool testMode = false;

/// Enabled further inter- and intra-solver synchronisation
bool syncMode = false;

namespace impl {

SolverInterfaceImpl:: SolverInterfaceImpl
(
  std::string participantName,
  int         accessorProcessRank,
  int         accessorCommunicatorSize,
  bool        serverMode )
:
  _accessorName(std::move(participantName)),
  _accessorProcessRank(accessorProcessRank),
  _accessorCommunicatorSize(accessorCommunicatorSize),
  _serverMode(serverMode)
{
  CHECK(!_accessorName.empty(), "Accessor has to be named!");
  CHECK(_accessorProcessRank >= 0, "Accessor process index has to be >= 0!");
  CHECK(_accessorCommunicatorSize >= 0, "Accessor process size has to be >= 0!");
  CHECK(_accessorProcessRank < _accessorCommunicatorSize,
        "Accessor process index has to be smaller than accessor process "
        << "size (given as " << _accessorProcessRank << ")!");

  /* When precice stops abruptly, e.g. an external solver crashes, the
     SolverInterfaceImpl destructor is never called. Since we still want
     to print the timings, we install the signal handler here. */
  // Disable SIGSEGV handler, because we don't want to interfere with crash backtrace.
  // signal(SIGSEGV, precice::utils::terminationSignalHandler);
  signal(SIGABRT, precice::utils::terminationSignalHandler);
  signal(SIGTERM, precice::utils::terminationSignalHandler);
  // SIGXCPU is emitted when the job is killed due to walltime limit on SuperMUC
  signal(SIGXCPU, precice::utils::terminationSignalHandler);
  // signal(SIGINT,  precice::utils::terminationSignalHandler);

  logging::setParticipant(_accessorName);
}

void SolverInterfaceImpl:: configure
(
  const std::string& configurationFileName )
{
  config::Configuration config;
  xml::configure(config.getXMLTag(), configurationFileName);
  if(_accessorProcessRank==0){
    INFO("This is preCICE version " << PRECICE_VERSION);
    INFO("Configuring preCICE with configuration: \"" << configurationFileName << "\"" );
  }
  configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceImpl:: configure
(
  const config::SolverInterfaceConfiguration& config )
{
  TRACE();

  Event e("configure"); // no precice::syncMode as this is not yet configured here
  utils::ScopedEventPrefix sep("configure/");

  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  Participant::resetParticipantCount();
  _meshLock.clear();

  _dimensions = config.getDimensions();
  _accessor = determineAccessingParticipant(config);

  CHECK(not (_accessor->useServer() && _accessor->useMaster()), "You cannot use a server and a master.");
  CHECK(_accessorCommunicatorSize==1 || _accessor->useMaster() || _accessor->useServer(),
        "A parallel participant needs either a master or a server communication configured");

  _clientMode = (not _serverMode) && _accessor->useServer();

  if(_accessor->useMaster()){
    utils::MasterSlave::configure(_accessorProcessRank, _accessorCommunicatorSize);
  }

  _participants = config.getParticipantConfiguration()->getParticipants();
  configureM2Ns(config.getM2NConfiguration());

  if (_serverMode){
    INFO("Run in server mode");
  }
  if (_clientMode){
    INFO("Run in client mode");
  }

  if (not _clientMode){
    INFO("Run in coupling mode");
    CHECK(_participants.size() > 1, "At least two participants need to be defined!");
    configurePartitions(config.getM2NConfiguration());
  }

  cplscheme::PtrCouplingSchemeConfiguration cplSchemeConfig =
      config.getCouplingSchemeConfiguration();
  _couplingScheme = cplSchemeConfig->getCouplingScheme(_accessorName);

  if (_serverMode || _clientMode){
    com::PtrCommunication com = _accessor->getClientServerCommunication();
    assertion(com.get() != nullptr);
    _requestManager = std::make_shared<RequestManager>(*this, com, _couplingScheme);
  }

  // Add meshIDs and data IDs
  for (const MeshContext* meshContext : _accessor->usedMeshContexts()) {
    const mesh::PtrMesh& mesh = meshContext->mesh;
    for (std::pair<std::string,int> nameID : mesh->getNameIDPairs()) {
      assertion(not utils::contained(nameID.first, _meshIDs));
      _meshIDs[nameID.first] = nameID.second;
    }
    assertion(_dataIDs.find(mesh->getID())==_dataIDs.end());
    _dataIDs[mesh->getID()] = std::map<std::string,int>();
    assertion(_dataIDs.find(mesh->getID())!=_dataIDs.end());
    for (const mesh::PtrData& data : mesh->data()) {
      assertion(_dataIDs[mesh->getID()].find(data->getName())==_dataIDs[mesh->getID()].end());
      _dataIDs[mesh->getID()][data->getName()] = data->getID();
    }
    std::string meshName = mesh->getName();
    mesh::PtrMeshConfiguration meshConfig = config.getMeshConfiguration();
  }
  // Register all MeshIds to the lock, but unlock them straight away as
  // writing is allowed after configuration.
  for (const auto& meshID : _meshIDs) {
      _meshLock.add(meshID.second, false);
  }
  
  utils::Parallel::initializeMPI(nullptr, nullptr);
  logging::setMPIRank(utils::Parallel::getProcessRank());
  utils::EventRegistry::instance().initialize("precice-" + _accessorName, "", utils::Parallel::getGlobalCommunicator());
  
  // Setup communication to server
  if (_clientMode){
    initializeClientServerCommunication();
  }
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()){
    initializeMasterSlaveCommunication();
  }

  auto & solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.start(precice::syncMode);
}

double SolverInterfaceImpl:: initialize()
{
  TRACE();
  auto & solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.pause(precice::syncMode);
  Event e("initialize", precice::syncMode);
  utils::ScopedEventPrefix sep("initialize/");
  
  if (_clientMode){
    DEBUG("Request perform initializations");
    _requestManager->requestInitialize();
  }
  else {
    // Setup communication

    INFO("Setting up master communication to coupling partner/s" );
    for (auto& m2nPair : _m2ns) {
        auto& bm2n = m2nPair.second;
        DEBUG((bm2n.isRequesting?"Awaiting master connection from ":"Establishing master connection to ") << bm2n.remoteName);
        bm2n.prepareEstablishment();
        bm2n.connectMasters();
        DEBUG("Established master connection " << (bm2n.isRequesting?"from ":"to ") << bm2n.remoteName);
    }
    INFO("Masters are connected");

    computePartitions();

    INFO("Setting up slaves communication to coupling partner/s" );
    for (auto& m2nPair : _m2ns) {
      auto& bm2n = m2nPair.second;
      DEBUG((bm2n.isRequesting?"Awaiting slaves connection from ":"Establishing slaves connection to ") << bm2n.remoteName);
      bm2n.connectSlaves();
      bm2n.cleanupEstablishment();
      DEBUG("Established slaves connection " << (bm2n.isRequesting?"from ":"to ") << bm2n.remoteName);
    }
    INFO("Slaves are connected");

    std::set<action::Action::Timing> timings;
    double dt = 0.0;

    for (PtrWatchPoint& watchPoint : _accessor->watchPoints()){
      watchPoint->initialize();
    }

    // Initialize coupling state, overwrite these values for restart
    double time = 0.0;
    int timestep = 1;

    _couplingScheme->initialize(time, timestep);

    dt = _couplingScheme->getNextTimestepMaxLength();

    timings.insert(action::Action::ALWAYS_POST);

    if (_couplingScheme->hasDataBeenExchanged()){
      timings.insert(action::Action::ON_EXCHANGE_POST);
      mapReadData();
    }

    performDataActions(timings, 0.0, 0.0, 0.0, dt);

    INFO(_couplingScheme->printCouplingState());
  }

  solverInitEvent.start(precice::syncMode);

  _meshLock.lockAll();

  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl:: initializeData ()
{
  TRACE();

  auto & solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.pause(precice::syncMode);

  Event e("initializeData", precice::syncMode);
  utils::ScopedEventPrefix sep("initializeData/");

  CHECK(_couplingScheme->isInitialized(),
        "initialize() has to be called before initializeData()");
  if (_clientMode){
    _requestManager->requestInitialzeData();
  }
  else {
    mapWrittenData();
    _couplingScheme->initializeData();
    double dt = _couplingScheme->getNextTimestepMaxLength();
    std::set<action::Action::Timing> timings;
    if (_couplingScheme->hasDataBeenExchanged()){
      timings.insert(action::Action::ON_EXCHANGE_POST);
      mapReadData();
    }
    performDataActions(timings, 0.0, 0.0, 0.0, dt);
    resetWrittenData();
    DEBUG("Plot output...");
    for (const io::ExportContext& context : _accessor->exportContexts()){
      if (context.timestepInterval != -1){
        std::ostringstream suffix;
        suffix << _accessorName << ".init";
        exportMesh(suffix.str());
        if (context.triggerSolverPlot){
          _couplingScheme->requireAction(std::string ("plot-output"));
        }
      }
    }
  }
  solverInitEvent.start(precice::syncMode);
}

double SolverInterfaceImpl:: advance
(
  double computedTimestepLength )
{
  TRACE(computedTimestepLength);

  // Events for the solver time, stopped when we enter, restarted when we leave advance
  auto & solverEvent = EventRegistry::instance().getStoredEvent("solver.advance");
  solverEvent.stop(precice::syncMode);
  auto & solverInitEvent = EventRegistry::instance().getStoredEvent("solver.initialize");
  solverInitEvent.stop(precice::syncMode);

  Event e("advance", precice::syncMode);
  utils::ScopedEventPrefix sep("advance/");

  CHECK(_couplingScheme->isInitialized(), "initialize() has to be called before advance()");
  _numberAdvanceCalls++;
  if (_clientMode){
    _requestManager->requestAdvance(computedTimestepLength);
  }
  else {
#   ifndef NDEBUG
    if(utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()){
      syncTimestep(computedTimestepLength);
    }
#   endif

    double timestepLength = 0.0; // Length of (full) current dt
    double timestepPart = 0.0;   // Length of computed part of (full) curr. dt
    double time = 0.0;


    // Update the coupling scheme time state. Necessary to get correct remainder.
    _couplingScheme->addComputedTime(computedTimestepLength);

    //double timestepLength = 0.0;
    if (_couplingScheme->hasTimestepLength()){
      timestepLength = _couplingScheme->getTimestepLength();
    }
    else {
      timestepLength = computedTimestepLength;
    }
    timestepPart = timestepLength - _couplingScheme->getThisTimestepRemainder();
    time = _couplingScheme->getTime();


    mapWrittenData();

    std::set<action::Action::Timing> timings;

    timings.insert(action::Action::ALWAYS_PRIOR);
    if (_couplingScheme->willDataBeExchanged(0.0)){
      timings.insert(action::Action::ON_EXCHANGE_PRIOR);
    }
    performDataActions(timings, time, computedTimestepLength, timestepPart, timestepLength);

    DEBUG("Advancing coupling scheme");
    _couplingScheme->advance();

    timings.clear();
    timings.insert(action::Action::ALWAYS_POST);
    if (_couplingScheme->hasDataBeenExchanged()){
      timings.insert(action::Action::ON_EXCHANGE_POST);
    }
    if (_couplingScheme->isCouplingTimestepComplete()){
      timings.insert(action::Action::ON_TIMESTEP_COMPLETE_POST);
    }
    performDataActions(timings, time, computedTimestepLength, timestepPart, timestepLength);

    if (_couplingScheme->hasDataBeenExchanged()){
      mapReadData();
    }

    INFO(_couplingScheme->printCouplingState());

    handleExports();

    // deactivated the reset of written data, as it deletes all data that is not communicated
    // within this cycle in the coupling data. This is not wanted forthe manifold mapping.
    //resetWrittenData();

  }
  _meshLock.lockAll();
  solverEvent.start(precice::syncMode);
  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl:: finalize()
{
  TRACE();

  // Events for the solver time, finally stopped here
  auto & solverEvent = EventRegistry::instance().getStoredEvent("solver.advance");
  solverEvent.stop(precice::syncMode);

  Event e("finalize"); // no precice::syncMode here as MPI is already finalized at destruction of this event
  utils::ScopedEventPrefix sep("finalize/");

  CHECK(_couplingScheme->isInitialized(), "initialize() has to be called before finalize()");
  _couplingScheme->finalize();
  _couplingScheme.reset();

  if (_clientMode){
    _requestManager->requestFinalize();
    _accessor->getClientServerCommunication()->closeConnection();
  }
  else {
    for (const io::ExportContext& context : _accessor->exportContexts()){
      if ( context.timestepInterval != -1 ){
        std::ostringstream suffix;
        suffix << _accessorName << ".final";
        exportMesh ( suffix.str() );
        if ( context.triggerSolverPlot ) {
          _couplingScheme->requireAction ( std::string ("plot-output") );
        }
      }
    }
    // Apply some final ping-pong to synch solver that run e.g. with a uni-directional coupling only
    // afterwards close connections
    std::string ping = "ping";
    std::string pong = "pong";
    for (auto &iter : _m2ns) {
      if( not utils::MasterSlave::isSlave()){
        if(iter.second.isRequesting){
          iter.second.m2n->getMasterCommunication()->send(ping,0);
          std::string receive = "init";
          iter.second.m2n->getMasterCommunication()->receive(receive,0);
          assertion(receive==pong);
        }
        else{
          std::string receive = "init";
          iter.second.m2n->getMasterCommunication()->receive(receive,0);
          assertion(receive==ping);
          iter.second.m2n->getMasterCommunication()->send(pong,0);
        }
      }
      iter.second.m2n->closeConnection();
    }
  }
  if(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()){
    utils::MasterSlave::_communication->closeConnection();
    utils::MasterSlave::_communication = nullptr;
  }

  if(_serverMode){
    _accessor->getClientServerCommunication()->closeConnection();
  }

  // Stop and print Event logging
  e.stop();
  utils::EventRegistry::instance().finalize();
  if (not precice::testMode and not precice::utils::MasterSlave::isSlave()) {
    utils::EventRegistry::instance().printAll();
  }

  // Tear down MPI and PETSc
  if (not precice::testMode && not _serverMode ) {
    utils::Petsc::finalize();
    utils::Parallel::finalizeMPI();
  }
  utils::Parallel::clearGroups();
  utils::EventRegistry::instance().clear();
}

int SolverInterfaceImpl:: getDimensions() const
{
  TRACE(_dimensions );
  return _dimensions;
}

bool SolverInterfaceImpl:: isCouplingOngoing() const
{
  TRACE();
  return _couplingScheme->isCouplingOngoing();
}

bool SolverInterfaceImpl:: isReadDataAvailable() const
{
  TRACE();
  return _couplingScheme->hasDataBeenExchanged();
}

bool SolverInterfaceImpl:: isWriteDataRequired
(
  double computedTimestepLength ) const
{
  TRACE(computedTimestepLength);
  return _couplingScheme->willDataBeExchanged(computedTimestepLength);
}

bool SolverInterfaceImpl:: isTimestepComplete() const
{
  TRACE();
  return _couplingScheme->isCouplingTimestepComplete();
}

bool SolverInterfaceImpl:: isActionRequired
(
  const std::string& action ) const
{
  TRACE(action, _couplingScheme->isActionRequired(action));
  return _couplingScheme->isActionRequired(action);
}

void SolverInterfaceImpl:: fulfilledAction
(
  const std::string& action )
{
  TRACE(action);
  if ( _clientMode ) {
    _requestManager->requestFulfilledAction(action);
  }
  _couplingScheme->performedAction(action);
}

bool SolverInterfaceImpl::hasToEvaluateSurrogateModel() const
{
 // std::cout<<"_isCoarseModelOptimizationActive() = "<<_couplingScheme->isCoarseModelOptimizationActive();
  return _couplingScheme->isCoarseModelOptimizationActive();
}

bool SolverInterfaceImpl::hasToEvaluateFineModel() const
{
  return not _couplingScheme->isCoarseModelOptimizationActive();
}

bool SolverInterfaceImpl:: hasMesh
(
  const std::string& meshName ) const
{
  TRACE(meshName);
  return utils::contained ( meshName, _meshIDs );
}

int SolverInterfaceImpl:: getMeshID
(
  const std::string& meshName ) const
{
  TRACE(meshName);
  const auto pos = _meshIDs.find(meshName);
  CHECK(pos != _meshIDs.end(), "Mesh with name \""<< meshName << "\" is not defined!" );
  return pos->second;
}

std::set<int> SolverInterfaceImpl:: getMeshIDs() const
{
  TRACE();
  std::set<int> ids;
  for (const impl::MeshContext* context : _accessor->usedMeshContexts()) {
    ids.insert ( context->mesh->getID() );
  }
  return ids;
}

bool SolverInterfaceImpl:: hasData
(
  const std::string& dataName, int meshID ) const
{
  TRACE(dataName, meshID );
  PRECICE_VALIDATE_MESH_ID(meshID);
  const auto & sub_dataIDs = _dataIDs.at(meshID);
  return sub_dataIDs.find(dataName)!= sub_dataIDs.end();
}

int SolverInterfaceImpl:: getDataID
(
  const std::string& dataName, int meshID ) const
{
  TRACE(dataName, meshID );
  PRECICE_VALIDATE_MESH_ID(meshID);
  CHECK(hasData(dataName, meshID),
        "Data with name \"" << dataName << "\" is not defined on mesh with ID \"" << meshID << "\".");
  return _dataIDs.at(meshID).at(dataName);
}

int SolverInterfaceImpl:: getMeshVertexSize
(
  int meshID ) const
{
  TRACE(meshID);
  int size = 0;
  if (_clientMode){
    size = _requestManager->requestGetMeshVertexSize(meshID);
  }
  else {
    PRECICE_REQUIRE_MESH_USE(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    assertion(context.mesh.get() != nullptr);
    size = context.mesh->vertices().size();
  }
  DEBUG("return " << size);
  return size;
}

/// @todo Currently not supported as we would need to re-compute the re-partition
void SolverInterfaceImpl:: resetMesh
(
  int meshID )
{
  TRACE(meshID);
  if (_clientMode){
    _requestManager->requestResetMesh(meshID);
  }
  else {
    PRECICE_VALIDATE_MESH_ID(meshID);
    impl::MeshContext& context = _accessor->meshContext(meshID);
    bool hasMapping = context.fromMappingContext.mapping
              || context.toMappingContext.mapping;
    bool isStationary =
          context.fromMappingContext.timing == mapping::MappingConfiguration::INITIAL &&
              context.toMappingContext.timing == mapping::MappingConfiguration::INITIAL;

    CHECK(!isStationary, "A mesh with only initial mappings  must not be reseted");
    CHECK(hasMapping, "A mesh with no mappings must not be reseted");

    DEBUG ( "Clear mesh positions for mesh \"" << context.mesh->getName() << "\"" );
    _meshLock.unlock(meshID);
    context.mesh->clear ();
  }
}

int SolverInterfaceImpl:: setMeshVertex
(
  int           meshID,
  const double* position )
{
  TRACE(meshID);
  Eigen::VectorXd internalPosition{
      Eigen::Map<const Eigen::VectorXd>{position, _dimensions}};
  DEBUG("Position = " << internalPosition);
  int index = -1;
  if ( _clientMode ){
    index = _requestManager->requestSetMeshVertex ( meshID, internalPosition);
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    DEBUG("MeshRequirement: " << context.meshRequirement);
    index = mesh->createVertex(internalPosition).getID();
    mesh->allocateDataValues();
  }
  return index;
}

void SolverInterfaceImpl:: setMeshVertices
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids )
{
  TRACE(meshID, size);
  if (_clientMode){
    _requestManager->requestSetMeshVertices(meshID, size, positions, ids);
  }
  else { //couplingMode
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    DEBUG("Set positions");
    const Eigen::Map<const Eigen::MatrixXd> posMatrix{
        positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
    for (int i=0; i < size; ++i){
      Eigen::VectorXd current(posMatrix.col(i));
      ids[i] = mesh->createVertex(current).getID();
    }
    mesh->allocateDataValues();
  }
}

void SolverInterfaceImpl:: getMeshVertices
(
  int        meshID,
  size_t     size,
  const int* ids,
  double*    positions ) const
{
  TRACE(meshID, size);
  if (_clientMode){
    _requestManager->requestGetMeshVertices(meshID, size, ids, positions);
  }
  else {
    PRECICE_REQUIRE_MESH_USE(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    DEBUG("Get positions");
    auto & vertices = mesh->vertices();
    assertion(size <= vertices.size(), size, vertices.size());
    Eigen::Map<Eigen::MatrixXd> posMatrix{
        positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
    for (size_t i=0; i < size; i++){
      const size_t id = ids[i];
      assertion(id < vertices.size(), id, vertices.size());
      posMatrix.col(i) = vertices[id].getCoords();
    }
  }
}

void SolverInterfaceImpl:: getMeshVertexIDsFromPositions (
  int           meshID,
  size_t        size,
  const double* positions,
  int*          ids ) const
{
  TRACE(meshID, size);
  if (_clientMode){
    _requestManager->requestGetMeshVertexIDsFromPositions(meshID, size, positions, ids);
  }
  else {
    PRECICE_REQUIRE_MESH_USE(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    DEBUG("Get IDs");
    const auto &vertices = mesh->vertices();
    assertion(vertices.size() <= size, vertices.size(), size);
    Eigen::Map<const Eigen::MatrixXd> posMatrix{
        positions, _dimensions, static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(size)};
    const auto vsize = vertices.size();
    for (size_t i = 0; i < size; i++) {
      size_t j=0;
      for (; j < vsize; j++) {
          if(math::equals(posMatrix.col(i), vertices[j].getCoords())) {
              break;
          }
      }
      CHECK(j != vsize, "Position " << i << "=" << ids[i] << " unknown!");
      ids[i] = j;
    }
  }
}



int SolverInterfaceImpl:: setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  TRACE(meshID, firstVertexID, secondVertexID );
  if ( _clientMode ){
    return _requestManager->requestSetMeshEdge ( meshID, firstVertexID, secondVertexID );
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::MeshRequirement::FULL ){
      DEBUG("Full mesh required.");
      mesh::PtrMesh& mesh = context.mesh;
      CHECK(mesh->isValidVertexID(firstVertexID),  " Given VertexID is invalid!");
      CHECK(mesh->isValidVertexID(secondVertexID), " Given VertexID is invalid!");
      mesh::Vertex& v0 = mesh->vertices()[firstVertexID];
      mesh::Vertex& v1 = mesh->vertices()[secondVertexID];
      return mesh->createEdge(v0, v1).getID ();
    }
  }
  return -1;
}


void SolverInterfaceImpl:: setMeshEdges
(
  int meshID,
  int size,
  const int* vertexIDs,
  int* edgeIDs )
{
  TRACE(meshID, size);
  if ( _clientMode ){
    _requestManager->requestSetMeshEdges( meshID, size, vertexIDs, edgeIDs);
    return;
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::MeshRequirement::FULL ){
      mesh::PtrMesh& mesh = context.mesh;
      CHECK(std::all_of(vertexIDs, vertexIDs + size * 2, [&mesh](int vid) {
                  return mesh->isValidVertexID(vid);
                  }),
              "Received invalid vertexIDs!");

      auto& vertices = mesh->vertices();
      for(int eidx = 0; eidx < size; ++eidx) {
          auto ida = vertexIDs[eidx*2];
          auto idb = vertexIDs[eidx*2+1];
          edgeIDs[eidx] = mesh->createEdge(vertices[ida], vertices[idb]).getID();
      }
    } else {
        std::fill(edgeIDs, edgeIDs+size, -1);
    }
  }
}

void SolverInterfaceImpl:: setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{
  TRACE(meshID, firstEdgeID,
                  secondEdgeID, thirdEdgeID );
  if ( _clientMode ){
    _requestManager->requestSetMeshTriangle ( meshID, firstEdgeID, secondEdgeID, thirdEdgeID );
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::MeshRequirement::FULL ){
      mesh::PtrMesh& mesh = context.mesh;
      CHECK(mesh->isValidEdgeID(firstEdgeID),  " Given EdgeID is invalid!");
      CHECK(mesh->isValidEdgeID(secondEdgeID), " Given EdgeID is invalid!");
      CHECK(mesh->isValidEdgeID(thirdEdgeID),  " Given EdgeID is invalid!");
      mesh::Edge& e0 = mesh->edges()[firstEdgeID];
      mesh::Edge& e1 = mesh->edges()[secondEdgeID];
      mesh::Edge& e2 = mesh->edges()[thirdEdgeID];
      mesh->createTriangle ( e0, e1, e2 );
    }
  }
}

void SolverInterfaceImpl:: setMeshTriangles
(
 int meshID,
 int size,
 const int* edgeIDs )
{
  TRACE(meshID, size);
  if ( _clientMode ){
    _requestManager->requestSetMeshTriangles ( meshID, size, edgeIDs );
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::MeshRequirement::FULL ){
      mesh::PtrMesh& mesh = context.mesh;
      CHECK(std::all_of( edgeIDs, edgeIDs+size*2, [&mesh](int vid) {
                  return mesh->isValidEdgeID(vid);
                  }), "Received invalid edgeIDs!");

      auto& edges = mesh->edges();
      for(int tidx = 0; tidx < size; ++tidx) {
          auto ida = edgeIDs[tidx*3];
          auto idb = edgeIDs[tidx*3+1];
          auto idc = edgeIDs[tidx*3+2];
          mesh->createTriangle(edges[ida], edges[idb], edges[idc]);
      }
    }
  }
}

void SolverInterfaceImpl:: setMeshTriangleWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID )
{
  TRACE(meshID, firstVertexID,
                secondVertexID, thirdVertexID);
  if (_clientMode){
    _requestManager->requestSetMeshTriangleWithEdges(meshID,
                                                     firstVertexID,
                                                     secondVertexID,
                                                     thirdVertexID);
    return;
  }
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext& context = _accessor->meshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL){
    mesh::PtrMesh& mesh = context.mesh;
    CHECK(mesh->isValidVertexID(firstVertexID),  " Given VertexID is invalid!");
    CHECK(mesh->isValidVertexID(secondVertexID), " Given VertexID is invalid!");
    CHECK(mesh->isValidVertexID(thirdVertexID),  " Given VertexID is invalid!");
    mesh::Vertex* vertices[3];
    vertices[0] = &mesh->vertices()[firstVertexID];
    vertices[1] = &mesh->vertices()[secondVertexID];
    vertices[2] = &mesh->vertices()[thirdVertexID];
    mesh::Edge* edges[3];
    edges[0] = & mesh->createUniqueEdge(*vertices[0], *vertices[1]);
    edges[1] = & mesh->createUniqueEdge(*vertices[1], *vertices[2]);
    edges[2] = & mesh->createUniqueEdge(*vertices[2], *vertices[0]);

    mesh->createTriangle(*edges[0], *edges[1], *edges[2]);
  }
}

void SolverInterfaceImpl:: setMeshQuad
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID,
  int fourthEdgeID )
{
  TRACE(meshID, firstEdgeID, secondEdgeID, thirdEdgeID,
                fourthEdgeID);
  if (_clientMode){
    _requestManager->requestSetMeshQuad(meshID, firstEdgeID, secondEdgeID,
                                        thirdEdgeID, fourthEdgeID);
  }
  else {
    PRECICE_REQUIRE_MESH_MODIFY(meshID);
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL){
      mesh::PtrMesh& mesh = context.mesh;
      CHECK(mesh->isValidEdgeID(firstEdgeID),  " Given EdgeID is invalid!");
      CHECK(mesh->isValidEdgeID(secondEdgeID), " Given EdgeID is invalid!");
      CHECK(mesh->isValidEdgeID(thirdEdgeID),  " Given EdgeID is invalid!");
      CHECK(mesh->isValidEdgeID(fourthEdgeID), " Given EdgeID is invalid!");
      mesh::Edge& e0 = mesh->edges()[firstEdgeID];
      mesh::Edge& e1 = mesh->edges()[secondEdgeID];
      mesh::Edge& e2 = mesh->edges()[thirdEdgeID];
      mesh::Edge& e3 = mesh->edges()[fourthEdgeID];
      mesh->createQuad(e0, e1, e2, e3);
    }
  }
}

void SolverInterfaceImpl:: setMeshQuadWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID,
  int fourthVertexID )
{
  TRACE(meshID, firstVertexID,
                secondVertexID, thirdVertexID, fourthVertexID);
  if (_clientMode){
    _requestManager->requestSetMeshQuadWithEdges(
        meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
    return;
  }
  PRECICE_REQUIRE_MESH_MODIFY(meshID);
  MeshContext& context = _accessor->meshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::MeshRequirement::FULL){
    mesh::PtrMesh& mesh = context.mesh;
    CHECK(mesh->isValidVertexID(firstVertexID),  " Given VertexID is invalid!");
    CHECK(mesh->isValidVertexID(secondVertexID), " Given VertexID is invalid!");
    CHECK(mesh->isValidVertexID(thirdVertexID),  " Given VertexID is invalid!");
    CHECK(mesh->isValidVertexID(fourthVertexID), " Given VertexID is invalid!");
    mesh::Vertex* vertices[4];
    vertices[0] = &mesh->vertices()[firstVertexID];
    vertices[1] = &mesh->vertices()[secondVertexID];
    vertices[2] = &mesh->vertices()[thirdVertexID];
    vertices[3] = &mesh->vertices()[fourthVertexID];
    mesh::Edge* edges[4];
    edges[0] = & mesh->createUniqueEdge(*vertices[0], *vertices[1]);
    edges[1] = & mesh->createUniqueEdge(*vertices[1], *vertices[2]);
    edges[2] = & mesh->createUniqueEdge(*vertices[2], *vertices[3]);
    edges[3] = & mesh->createUniqueEdge(*vertices[3], *vertices[0]);

    mesh->createQuad(*edges[0], *edges[1], *edges[2], *edges[3]);
  }
}

void SolverInterfaceImpl:: mapWriteDataFrom
(
  int fromMeshID )
{
  TRACE(fromMeshID);
  if (_clientMode){
    _requestManager->requestMapWriteDataFrom(fromMeshID);
    return;
  }
  PRECICE_VALIDATE_MESH_ID(fromMeshID);
  impl::MeshContext& context = _accessor->meshContext(fromMeshID);
  impl::MappingContext& mappingContext = context.fromMappingContext;
  if (mappingContext.mapping.use_count() == 0){
    ERROR("From mesh \"" << context.mesh->getName()
                   << "\", there is no mapping defined");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    DEBUG("Compute mapping from mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  for (impl::DataContext& context : _accessor->writeDataContexts()) {
    if (context.mesh->getID() == fromMeshID){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      DEBUG("Map data \"" << context.fromData->getName()
                   << "\" from mesh \"" << context.mesh->getName() << "\"");
      assertion(mappingContext.mapping==context.mappingContext.mapping);
      mappingContext.mapping->map(inDataID, outDataID);
    }
  }
  mappingContext.hasMappedData = true;
}


void SolverInterfaceImpl:: mapReadDataTo
(
  int toMeshID )
{
  TRACE(toMeshID);
  if (_clientMode){
    _requestManager->requestMapReadDataTo(toMeshID);
    return;
  }
  PRECICE_VALIDATE_MESH_ID(toMeshID);
  impl::MeshContext& context = _accessor->meshContext(toMeshID);
  impl::MappingContext& mappingContext = context.toMappingContext;
  if (mappingContext.mapping.use_count() == 0){
    ERROR("From mesh \"" << context.mesh->getName()
                   << "\", there is no mapping defined!");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    DEBUG("Compute mapping from mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  for (impl::DataContext& context : _accessor->readDataContexts()) {
    if (context.mesh->getID() == toMeshID){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      DEBUG("Map data \"" << context.fromData->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      assertion(mappingContext.mapping==context.mappingContext.mapping);
      mappingContext.mapping->map(inDataID, outDataID);
      DEBUG("First mapped values = " << utils::firstN(context.toData->values(), 10));
    }
  }
  mappingContext.hasMappedData = true;
}

void SolverInterfaceImpl:: writeBlockVectorData
(
  int           fromDataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  TRACE(fromDataID, size);
  PRECICE_VALIDATE_DATA_ID(fromDataID);
  if (size == 0)
    return;
  assertion(valueIndices != nullptr);
  assertion(values != nullptr);
  if (_clientMode){
    _requestManager->requestWriteBlockVectorData(fromDataID, size, valueIndices, values);
  }
  else { //couplingMode
    PRECICE_REQUIRE_DATA_WRITE(fromDataID);
    DataContext& context = _accessor->dataContext(fromDataID);
    CHECK(context.fromData->getDimensions()==_dimensions,
        "You cannot call writeBlockVectorData on the scalar data type " << context.fromData->getName());
    assertion(context.toData.get() != nullptr);
    auto& valuesInternal = context.fromData->values();
    for (int i=0; i < size; i++){
      const auto valueIndex = valueIndices[i];
      CHECK(0 <= valueIndex && valueIndex < valuesInternal.size()/context.fromData->getDimensions(), "Value index out of range");
      int offsetInternal = valueIndex*_dimensions;
      int offset = i*_dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        assertion(offset+dim < valuesInternal.size(),
                   offset+dim, valuesInternal.size());
        valuesInternal[offsetInternal + dim] = values[offset + dim];
      }
    }
  }
}

void SolverInterfaceImpl:: writeVectorData
(
  int           fromDataID,
  int           valueIndex,
  const double* value )
{
  TRACE(fromDataID, valueIndex );
  PRECICE_VALIDATE_DATA_ID(fromDataID);
  DEBUG("value = " << Eigen::Map<const Eigen::VectorXd>(value, _dimensions));
  if (_clientMode){
    _requestManager->requestWriteVectorData(fromDataID, valueIndex, value);
  }
  else {
    PRECICE_REQUIRE_DATA_WRITE(fromDataID);
    DataContext& context = _accessor->dataContext(fromDataID);
    CHECK(context.fromData->getDimensions()==_dimensions,
        "You cannot call writeVectorData on the scalar data type " << context.fromData->getName());
    assertion(context.toData.get() != nullptr);
    auto& values = context.fromData->values();
    CHECK(0 <= valueIndex && valueIndex < values.size()/context.fromData->getDimensions(), "Value index out of range");
    int offset = valueIndex * _dimensions;
    for (int dim=0; dim < _dimensions; dim++){
      values[offset+dim] = value[dim];
    }

  }
}

void SolverInterfaceImpl:: writeBlockScalarData
(
  int           fromDataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  TRACE(fromDataID, size);
  PRECICE_VALIDATE_DATA_ID(fromDataID);
  if (size == 0)
    return;
  assertion(valueIndices != nullptr);
  assertion(values != nullptr);
  if (_clientMode){
    _requestManager->requestWriteBlockScalarData(fromDataID, size, valueIndices, values);
  }
  else {
    PRECICE_REQUIRE_DATA_WRITE(fromDataID);
    DataContext& context = _accessor->dataContext(fromDataID);
    CHECK(context.fromData->getDimensions()==1,
        "You cannot call writeBlockScalarData on the vector data type " << context.fromData->getName());
    assertion(context.toData.get() != nullptr);
    auto& valuesInternal = context.fromData->values();
    for (int i=0; i < size; i++){
      const auto valueIndex = valueIndices[i];
      CHECK(0 <= valueIndex && valueIndex < valuesInternal.size()/context.fromData->getDimensions(), "Value index out of range");
      assertion(i < valuesInternal.size(), i, valuesInternal.size());
      valuesInternal[valueIndex] = values[i];
    }
  }
}

void SolverInterfaceImpl:: writeScalarData
(
  int    fromDataID,
  int    valueIndex,
  double value)
{
  TRACE(fromDataID, valueIndex, value );
  PRECICE_VALIDATE_DATA_ID(fromDataID);
  CHECK(valueIndex >= -1, "Invalid value index (" << valueIndex << ") when writing scalar data!");
  if (_clientMode){
    _requestManager->requestWriteScalarData(fromDataID, valueIndex, value);
  }
  else {
    PRECICE_REQUIRE_DATA_WRITE(fromDataID);
    DataContext& context = _accessor->dataContext(fromDataID);
    CHECK(context.fromData->getDimensions()==1,
        "You cannot call writeScalarData on the vector data type " << context.fromData->getName());
    assertion(context.toData);
    auto& values = context.fromData->values();
    CHECK(0 <= valueIndex && valueIndex < values.size()/context.fromData->getDimensions(), "Value index out of range");
    values[valueIndex] = value;
  }
}

void SolverInterfaceImpl:: readBlockVectorData
(
  int        toDataID,
  int        size,
  const int* valueIndices,
  double*    values ) const
{
  TRACE(toDataID, size);
  PRECICE_VALIDATE_DATA_ID(toDataID);
  if (size == 0)
    return;
  assertion(valueIndices != nullptr);
  assertion(values != nullptr);
  if (_clientMode){
    _requestManager->requestReadBlockVectorData(toDataID, size, valueIndices, values);
  }
  else { //couplingMode
    PRECICE_REQUIRE_DATA_READ(toDataID);
    DataContext& context = _accessor->dataContext(toDataID);
    CHECK(context.toData->getDimensions()==_dimensions,
        "You cannot call readBlockVectorData on the scalar data type " << context.toData->getName());
    assertion(context.fromData.get() != nullptr);
    auto& valuesInternal = context.toData->values();
    for (int i=0; i < size; i++){
      const auto valueIndex = valueIndices[i];
      CHECK(0 <= valueIndex && valueIndex < valuesInternal.size()/context.fromData->getDimensions(), "Value index out of range");
      int offsetInternal = valueIndex * _dimensions;
      int offset = i * _dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        assertion(offsetInternal+dim < valuesInternal.size(),
                   offsetInternal+dim, valuesInternal.size());
        values[offset + dim] = valuesInternal[offsetInternal + dim];
      }
    }
  }
}

void SolverInterfaceImpl:: readVectorData
(
  int     toDataID,
  int     valueIndex,
  double* value ) const
{
  TRACE(toDataID, valueIndex);
  PRECICE_VALIDATE_DATA_ID(toDataID);
  CHECK(valueIndex >= -1, "Invalid value index ( " << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadVectorData(toDataID, valueIndex, value);
  }
  else {
    PRECICE_REQUIRE_DATA_READ(toDataID);
    DataContext& context = _accessor->dataContext(toDataID);
    CHECK(context.toData->getDimensions()==_dimensions,
        "You cannot call readVectorData on the scalar data type " << context.toData->getName());
    assertion(context.fromData);
    auto& values = context.toData->values();
    CHECK(0 <= valueIndex && valueIndex < values.size()/context.fromData->getDimensions(), "Value index out of range");
    int offset = valueIndex * _dimensions;
    for (int dim=0; dim < _dimensions; dim++){
      value[dim] = values[offset + dim];
    }

  }
  DEBUG("read value = " << Eigen::Map<const Eigen::VectorXd>(value, _dimensions));
}

void SolverInterfaceImpl:: readBlockScalarData
(
  int        toDataID,
  int        size,
  const int* valueIndices,
  double*    values ) const
{
  TRACE(toDataID, size);
  PRECICE_VALIDATE_DATA_ID(toDataID);
  if (size == 0)
    return;
  DEBUG("size = " << size);
  assertion(valueIndices != nullptr);
  assertion(values != nullptr);
  if (_clientMode){
    _requestManager->requestReadBlockScalarData(toDataID, size, valueIndices, values);
  }
  else {
    PRECICE_REQUIRE_DATA_READ(toDataID);
    DataContext& context = _accessor->dataContext(toDataID);
    CHECK(context.toData->getDimensions()==1,
        "You cannot call readBlockScalarData on the vector data type " << context.toData->getName());
    assertion(context.fromData.get() != nullptr);
    auto& valuesInternal = context.toData->values();
    for (int i=0; i < size; i++){
        const auto valueIndex = valueIndices[i];
        CHECK(0 <= valueIndex && valueIndex < valuesInternal.size(), "Value index out of range");
        values[i] = valuesInternal[valueIndex];
    }
  }
}

void SolverInterfaceImpl:: readScalarData
(
  int     toDataID,
  int     valueIndex,
  double& value ) const
{
  TRACE(toDataID, valueIndex, value);
  PRECICE_VALIDATE_DATA_ID(toDataID);
  CHECK(valueIndex >= -1, "Invalid value index ( " << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadScalarData(toDataID, valueIndex, value);
  }
  else {
    PRECICE_REQUIRE_DATA_READ(toDataID);
    DataContext& context = _accessor->dataContext(toDataID);
    CHECK(context.toData->getDimensions()==1,
        "You cannot call readScalarData on the vector data type " << context.toData->getName());
    assertion(context.fromData);
    auto& values = context.toData->values();
    CHECK(0 <= valueIndex && valueIndex < values.size(), "Value index out of range");
    value = values[valueIndex];

  }
  DEBUG("Read value = " << value);
}

void SolverInterfaceImpl:: exportMesh
(
  const std::string& filenameSuffix,
  int                exportType ) const
{
  TRACE(filenameSuffix, exportType );
  // Export meshes
  //const ExportContext& context = _accessor->exportContext();
  for (const io::ExportContext& context : _accessor->exportContexts()) {
    DEBUG ( "Export type = " << exportType );
    bool exportAll = exportType == io::constants::exportAll();
    bool exportThis = context.exporter->getType() == exportType;
    if ( exportAll || exportThis ){
      for (const MeshContext* meshContext : _accessor->usedMeshContexts()) {
        std::string name = meshContext->mesh->getName() + "-" + filenameSuffix;
        DEBUG ( "Exporting mesh to file \"" << name << "\" at location \"" << context.location << "\"" );
        context.exporter->doExport ( name, context.location, *(meshContext->mesh) );
      }
    }
  }
}


MeshHandle SolverInterfaceImpl:: getMeshHandle
(
  const std::string& meshName )
{
  TRACE(meshName);
  assertion(not _clientMode);
  for (MeshContext* context : _accessor->usedMeshContexts()){
    if (context->mesh->getName() == meshName){
      return {context->mesh->content()};
    }
  }
  ERROR("Participant \"" << _accessorName
               << "\" does not use mesh \"" << meshName << "\"!");
}

void SolverInterfaceImpl:: runServer()
{
  assertion(_serverMode);
  initializeClientServerCommunication();
  _requestManager->handleRequests();
}

void SolverInterfaceImpl:: configureM2Ns
(
  const m2n::M2NConfiguration::SharedPointer& config )
{
  TRACE();
  for (const auto& m2nTuple : config->m2ns()) {
    std::string comPartner("");
    bool isRequesting = false;
    if (std::get<1>(m2nTuple) == _accessorName){
      comPartner = std::get<2>(m2nTuple);
      isRequesting = true;
    }
    else if (std::get<2>(m2nTuple) == _accessorName){
      comPartner = std::get<1>(m2nTuple);
    }
    if (not comPartner.empty()){
      for (const impl::PtrParticipant& participant : _participants) {
        if (participant->getName() == comPartner){
          if (participant->useServer()){
            comPartner += "Server";
          }
          assertion(not utils::contained(comPartner, _m2ns), comPartner);
          assertion(std::get<0>(m2nTuple));

          _m2ns[comPartner] = [&]{
              m2n::BoundM2N bound;
              bound.m2n = std::get<0>(m2nTuple);
              bound.localName = _accessorName;
              bound.remoteName = comPartner;
              bound.isRequesting = isRequesting;
              bound.localServer = _serverMode;
              return bound;
          }();
        }
      }
    }
  }
}

void SolverInterfaceImpl:: configurePartitions
(
  const m2n::M2NConfiguration::SharedPointer& m2nConfig )
{
  TRACE();
  for (MeshContext* context : _accessor->usedMeshContexts()) {
    if ( context->provideMesh ) { // Accessor provides mesh
      CHECK ( context->receiveMeshFrom.empty(),
              "Participant \"" << _accessorName << "\" cannot provide "
              << "and receive mesh " << context->mesh->getName() << "!" );

      context->partition = partition::PtrPartition(new partition::ProvidedPartition(context->mesh));

      for (auto& receiver : _participants ) {
        for (auto& receiverContext : receiver->usedMeshContexts()) {
          if(receiverContext->receiveMeshFrom == _accessorName && receiverContext->mesh->getName() == context->mesh->getName()){
            //CHECK( not hasToSend, "Mesh " << context->mesh->getName() << " can currently only be received once.")

            // meshRequirement has to be copied from "from" to provide", since
            // mapping are only defined at "provide"
            if(receiverContext->meshRequirement > context->meshRequirement){
              context->meshRequirement = receiverContext->meshRequirement;
            }
            m2n::PtrM2N m2n = m2nConfig->getM2N( receiver->getName(), _accessorName );
            m2n->createDistributedCommunication(context->mesh);
            context->partition->addM2N(m2n);
          }
        }
      }
      /// @todo support offset??

    }
    else { // Accessor receives mesh
      CHECK(not context->receiveMeshFrom.empty(),
            "Participant \"" << _accessorName << "\" must either provide or receive the mesh " << context->mesh->getName() << "!")
      CHECK(not context->provideMesh,
            "Participant \"" << _accessorName << "\" cannot provide and receive mesh " << context->mesh->getName() << "!" );
      std::string receiver ( _accessorName );
      std::string provider ( context->receiveMeshFrom );
      DEBUG ( "Receiving mesh from " << provider );
      
      context->partition = partition::PtrPartition(new partition::ReceivedPartition(context->mesh, context->geoFilter, context->safetyFactor));

      m2n::PtrM2N m2n = m2nConfig->getM2N ( receiver, provider );
      m2n->createDistributedCommunication(context->mesh);
      context->partition->addM2N(m2n);
      context->partition->setFromMapping(context->fromMappingContext.mapping);
      context->partition->setToMapping(context->toMappingContext.mapping);
    }
  }
}

void SolverInterfaceImpl:: computePartitions()
{
  //We need to do this in two loops: First, communicate the mesh and later compute the partition.
  //Originally, this was done in one loop. This however gave deadlock if two meshes needed to be communicated cross-wise.
  //Both loops need a different sorting
  auto &contexts = _accessor->usedMeshContexts();

  // sort meshContexts by name, for communication in right order.
  std::sort(contexts.begin(), contexts.end(),
            [](MeshContext const *const lhs, MeshContext const *const rhs) -> bool {
              return lhs->mesh->getName() < rhs->mesh->getName();
            });

  for (MeshContext *meshContext : contexts) {
    meshContext->partition->communicate();
  }

  // pull provided meshes up front, to have them ready for the decomposition
  std::stable_partition(contexts.begin(), contexts.end(),
                        [](MeshContext const *const meshContext) -> bool {
                          return meshContext->provideMesh;
                        });

  for (MeshContext *meshContext : contexts) {
    meshContext->partition->compute();
    meshContext->mesh->computeState();
    meshContext->mesh->allocateDataValues();
  }
}


void SolverInterfaceImpl:: mapWrittenData()
{
  TRACE();
  using namespace mapping;
  MappingConfiguration::Timing timing;
  // Compute mappings
  for (impl::MappingContext& context : _accessor->writeMappingContexts()) {
    timing = context.timing;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    bool hasComputed = context.mapping->hasComputedMapping();
    if (rightTime && not hasComputed){
      INFO("Compute write mapping from mesh \""
          << _accessor->meshContext(context.fromMeshID).mesh->getName()
          << "\" to mesh \""
          << _accessor->meshContext(context.toMeshID).mesh->getName()
          << "\".");

      context.mapping->computeMapping();
    }
  }

  // Map data
  for (impl::DataContext& context : _accessor->writeDataContexts()) {
    timing = context.mappingContext.timing;
    bool hasMapping = context.mappingContext.mapping.get() != nullptr;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    bool hasMapped = context.mappingContext.hasMappedData;
    if (hasMapping && rightTime && (not hasMapped)){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      DEBUG("Map data \"" << context.fromData->getName()
                   << "\" from mesh \"" << context.mesh->getName() << "\"");
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      DEBUG("Map from dataID " << inDataID << " to dataID: " << outDataID);
      context.mappingContext.mapping->map(inDataID, outDataID);
      DEBUG("First mapped values = " << utils::firstN(context.toData->values(), 10));
    }
  }

  // Clear non-stationary, non-incremental mappings
  for (impl::MappingContext& context : _accessor->writeMappingContexts()) {
    bool isStationary = context.timing
                        == MappingConfiguration::INITIAL;
    if (not isStationary){
        context.mapping->clear();
    }
    context.hasMappedData = false;
  }
}

void SolverInterfaceImpl:: mapReadData()
{
  TRACE();
  mapping::MappingConfiguration::Timing timing;
  // Compute mappings
  for (impl::MappingContext& context : _accessor->readMappingContexts()) {
    timing = context.timing;
    bool mapNow = timing == mapping::MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == mapping::MappingConfiguration::INITIAL;
    bool hasComputed = context.mapping->hasComputedMapping();
    if (mapNow && not hasComputed){
      INFO("Compute read mapping from mesh \""
              << _accessor->meshContext(context.fromMeshID).mesh->getName()
              << "\" to mesh \""
              << _accessor->meshContext(context.toMeshID).mesh->getName()
              << "\".");

      context.mapping->computeMapping();
    }
  }

  // Map data
  for (impl::DataContext& context : _accessor->readDataContexts()) {
    timing = context.mappingContext.timing;
    bool mapNow = timing == mapping::MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == mapping::MappingConfiguration::INITIAL;
    bool hasMapping = context.mappingContext.mapping.get() != nullptr;
    bool hasMapped = context.mappingContext.hasMappedData;
    if (mapNow && hasMapping && (not hasMapped)){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
      DEBUG("Map read data \"" << context.fromData->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      context.mappingContext.mapping->map(inDataID, outDataID);
      DEBUG("First mapped values = " << utils::firstN(context.toData->values(), 10));
    }
  }

  // Clear non-initial, non-incremental mappings
  for (impl::MappingContext& context : _accessor->readMappingContexts()) {
    bool isStationary = context.timing
              == mapping::MappingConfiguration::INITIAL;
    if (not isStationary){
      context.mapping->clear();
    }
    context.hasMappedData = false;
  }
}

void SolverInterfaceImpl:: performDataActions
(
  const std::set<action::Action::Timing>& timings,
  double                 time,
  double                 dt,
  double                 partFullDt,
  double                 fullDt )
{
  TRACE();
  assertion(not _clientMode);
  for (action::PtrAction& action : _accessor->actions()) {
    if (timings.find(action->getTiming()) != timings.end()){
      action->performAction(time, dt, partFullDt, fullDt);
    }
  }
}

void SolverInterfaceImpl:: handleExports()
{
  TRACE();
  assertion(not _clientMode);
  //timesteps was already incremented before
  int timesteps = _couplingScheme->getTimesteps()-1;

  for (const io::ExportContext& context : _accessor->exportContexts()) {
    if (_couplingScheme->isCouplingTimestepComplete() || context.everyIteration){
      if (context.timestepInterval != -1){
        if (timesteps % context.timestepInterval == 0){
          if (context.everyIteration){
            std::ostringstream everySuffix;
            everySuffix << _accessorName << ".it" << _numberAdvanceCalls;
            exportMesh(everySuffix.str());
          }
          std::ostringstream suffix;
          suffix << _accessorName << ".dt" << _couplingScheme->getTimesteps()-1;
          exportMesh(suffix.str());
          if (context.triggerSolverPlot){
            _couplingScheme->requireAction( std::string ("plot-output") );
          }
        }
      }
    }
  }

  if (_couplingScheme->isCouplingTimestepComplete()){
    // Export watch point data
    for (const PtrWatchPoint& watchPoint : _accessor->watchPoints()) {
      watchPoint->exportPointData(_couplingScheme->getTime());
    }
  }
}

void SolverInterfaceImpl:: resetWrittenData()
{
  TRACE();
  for (DataContext& context : _accessor->writeDataContexts()) {
    context.fromData->values() = Eigen::VectorXd::Zero(context.fromData->values().size());
    if (context.toData != context.fromData){
      context.toData->values() = Eigen::VectorXd::Zero(context.toData->values().size());
    }
  }
}

PtrParticipant SolverInterfaceImpl:: determineAccessingParticipant
(
   const config::SolverInterfaceConfiguration& config )
{
  const auto& partConfig = config.getParticipantConfiguration();
  for (const PtrParticipant& participant : partConfig->getParticipants()) {
    if ( participant->getName() == _accessorName ) {
      return participant;
    }
  }
  ERROR("Accessing participant \"" << _accessorName << "\" is not defined in configuration!");
}


void SolverInterfaceImpl:: initializeClientServerCommunication()
{
  TRACE();
  com::PtrCommunication com = _accessor->getClientServerCommunication();
  assertion(com.get() != nullptr);
  if ( _serverMode ){
    INFO("Setting up communication to client" );
    com->acceptConnection ( _accessorName + "Server", _accessorName, 0);
  }
  else {
    INFO("Setting up communication to server" );
    com->requestConnection( _accessorName + "Server", _accessorName,
                            _accessorProcessRank, _accessorCommunicatorSize );
  }
}

void SolverInterfaceImpl:: initializeMasterSlaveCommunication()
{
  TRACE();

  Event e("com.initializeMasterSlaveCom", precice::syncMode);
  //slaves create new communicator with ranks 0 to size-2
  //therefore, the master uses a rankOffset and the slaves have to call request
  // with that offset
  int rankOffset = 1;
  if ( utils::MasterSlave::isMaster() ){
    INFO("Setting up communication to slaves" );
    utils::MasterSlave::_communication->acceptConnection ( _accessorName + "Master", _accessorName, utils::MasterSlave::getRank());
    utils::MasterSlave::_communication->setRankOffset(rankOffset);
  }
  else {
    assertion(utils::MasterSlave::isSlave());
    utils::MasterSlave::_communication->requestConnection( _accessorName + "Master", _accessorName,
                            _accessorProcessRank-rankOffset, _accessorCommunicatorSize-rankOffset );
  }
}

void SolverInterfaceImpl:: syncTimestep(double computedTimestepLength)
{
  assertion(utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave());
  if(utils::MasterSlave::isSlave()){
    utils::MasterSlave::_communication->send(computedTimestepLength, 0);
  }
  else if(utils::MasterSlave::isMaster()){
    for(int rankSlave = 1; rankSlave < _accessorCommunicatorSize; rankSlave++){
      double dt;
      utils::MasterSlave::_communication->receive(dt, rankSlave);
      CHECK(math::equals(dt, computedTimestepLength),
            "Ambiguous timestep length when calling request advance from several processes!");
    }
  }
}

}} // namespace precice, impl

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverInterfaceImpl.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "precice/impl/RequestManager.hpp"
#include "precice/config/Configuration.hpp"
#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "tarch/la/WrappedVector.h"
#include "tarch/logging/CommandLineLogger.h"
#include "precice/config/LogFilterConfiguration.hpp"
#include "precice/config/LogOutputFormatConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Group.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Merge.hpp"
#include "io/ExportVTK.hpp"
#include "io/ExportVRML.hpp"
#include "io/ExportContext.hpp"
#include "io/SimulationStateIO.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "query/FindClosest.hpp"
#include "query/FindVoxelContent.hpp"
#include "spacetree/config/SpacetreeConfiguration.hpp"
#include "spacetree/Spacetree.hpp"
#include "spacetree/ExportSpacetree.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/Constants.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointCommunication.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "geometry/Geometry.hpp"
#include "geometry/ImportGeometry.hpp"
#include "geometry/CommunicatedGeometry.hpp"
#include "geometry/SolverGeometry.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "utils/MasterSlave.hpp"
#include "mapping/Mapping.hpp"
#include <set>
#include <limits>
#include <cstring>
#include <algorithm>
#include "utils/EventTimings.hpp"
#include "boost/tuple/tuple.hpp"

#include <signal.h> // used for installing crash handler

#ifndef PRECICE_NO_PETSC
#include "petsc.h"
#endif


using precice::utils::Event;

namespace precice {

bool testMode = false;

namespace impl {

tarch::logging::Log SolverInterfaceImpl::
  _log ("precice::impl::SolverInterfaceImpl");

SolverInterfaceImpl:: SolverInterfaceImpl
(
  const std::string& accessorName,
  int                accessorProcessRank,
  int                accessorCommunicatorSize,
  bool               serverMode )
:
  _accessorName(accessorName),
  _accessorProcessRank(accessorProcessRank),
  _accessorCommunicatorSize(accessorCommunicatorSize),
  _accessor(),
  _dimensions(0),
  _geometryMode(false),
  _restartMode(false),
  _serverMode(serverMode),
  _clientMode(false),
  _meshIDs(),
  _dataIDs(),
  _exportVTKNeighbors(),
  _m2ns(),
  _participants(),
  _checkpointTimestepInterval(-1),
  _checkpointFileName("precice_checkpoint_" + _accessorName),
  _numberAdvanceCalls(0),
  _requestManager(NULL)
{
  preciceCheck(_accessorProcessRank >= 0, "SolverInterfaceImpl()",
               "Accessor process index has to be >= 0!");
  preciceCheck(_accessorCommunicatorSize >= 0, "SolverInterfaceImpl()",
               "Accessor process size has to be >= 0!");
  preciceCheck(_accessorProcessRank < _accessorCommunicatorSize,
               "SolverInterfaceImpl()",
               "Accessor process index has to be smaller than accessor process "
               << "size (given as " << _accessorProcessRank << ")!");
  precice::utils::Events_Init();

  /* When precice stops abruptly, e.g. an external solver crashes, the
     SolverInterfaceImpl destructor is never called. Since we still want
     to print the timings, we install the signal handler here. */
  signal(SIGSEGV, precice::utils::EventRegistry::signal_handler);
  signal(SIGABRT, precice::utils::EventRegistry::signal_handler);
  signal(SIGTERM, precice::utils::EventRegistry::signal_handler);
}

SolverInterfaceImpl:: ~SolverInterfaceImpl()
{
  if (_requestManager != NULL){
    delete _requestManager;
  }
  precice::utils::Events_Finalize();
  if (not precice::utils::MasterSlave::_slaveMode) {
    precice::utils::EventRegistry r;
    r.print();
    r.print("EventTimings.log", true);
  }
}

void SolverInterfaceImpl:: configure
(
  const std::string& configurationFileName )
{
  typedef tarch::logging::CommandLineLogger Logger;
  // By default, debugging is turned on with a filter list entry. This removes
  // entry and turns off all debug messages until configuration.
  Logger::getInstance().clearFilterList();
  Logger::FilterListEntry filter("", true); // All off
  Logger::getInstance().addFilterListEntry(filter);

  preciceTrace1("configure()", configurationFileName );
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  Participant::resetParticipantCount();

  config::Configuration config;
  utils::configure(config.getXMLTag(), configurationFileName);
  //preciceCheck ( config.isValid(), "configure()", "Invalid configuration file!" );

  const config::LogFilterConfiguration& logFilterConfig =
      config.getLogFilterConfiguration();
  Logger::getInstance().clearFilterList();
  Logger::getInstance().addFilterListEntries(logFilterConfig.getFilterList());

  const config::LogOutputFormatConfiguration& logFormatConfig =
      config.getLogFormatConfiguration();
  Logger::getInstance().setLogFormat(
      logFormatConfig.getLogColumnSeparator(),
      logFormatConfig.getLogTimeStamp(),
      logFormatConfig.getLogTimeStampHumanReadable(),
      logFormatConfig.getLogMachineName(),
      logFormatConfig.getLogMessageType(),
      logFormatConfig.getLogTrace(),
      "");
  configure(config.getSolverInterfaceConfiguration());
}

void SolverInterfaceImpl:: configure
(
  const config::SolverInterfaceConfiguration& config )
{
  preciceTrace("configure()");
  _dimensions = config.getDimensions();
  _geometryMode = config.isGeometryMode ();
  _restartMode = config.isRestartMode ();
  _accessor = determineAccessingParticipant(config);

  preciceCheck(not (_accessor->useServer() && _accessor->useMaster()),
                     "configure()", "You cannot use a server and a master.");
  preciceCheck(not (_restartMode && _accessor->useMaster()),"configure()",
                      "To restart while using a master is not yet supported");

  _clientMode = (not _serverMode) && _accessor->useServer();

  if(_accessor->useMaster()){
    utils::MasterSlave::configure(_accessorProcessRank, _accessorCommunicatorSize);
  }

  _participants = config.getParticipantConfiguration()->getParticipants();
  configureM2Ns(config.getM2NConfiguration());

  if (_serverMode){
    preciceInfo("configure()", "[PRECICE] Run in server mode");
  }
  if (_clientMode){
    preciceInfo("configure()", "[PRECICE] Run in client mode");
  }

  if (_geometryMode){
    preciceInfo("configure()", "[PRECICE] Run in geometry mode");
    preciceCheck(_participants.size() == 1, "configure()",
                 "Only one participant can be defined in geometry mode!");
    configureSolverGeometries(config.getM2NConfiguration());
  }
  else if (not _clientMode){
    preciceInfo("configure()", "[PRECICE] Run in coupling mode");
    preciceCheck(_participants.size() > 1,
                 "configure()", "At least two participants need to be defined!");
    configureSolverGeometries(config.getM2NConfiguration());
  }

  // Set coupling scheme. In geometry mode, an uncoupled scheme is automatically
  // created.
  cplscheme::PtrCouplingSchemeConfiguration cplSchemeConfig =
      config.getCouplingSchemeConfiguration();
  _couplingScheme = cplSchemeConfig->getCouplingScheme(_accessorName);

  if (_restartMode){
    _couplingScheme->requireAction(constants::actionReadSimulationCheckpoint());
  }

  if (_serverMode || _clientMode){
    com::Communication::SharedPointer com = _accessor->getClientServerCommunication();
    assertion(com.get() != NULL);
    _requestManager = new RequestManager(_geometryMode, *this, com, _couplingScheme);
  }

  // Add meshIDs, data IDs, and spacetrees
  for (MeshContext* meshContext : _accessor->usedMeshContexts()) {
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
    spacetree::PtrSpacetreeConfiguration spacetreeConfig = config.getSpacetreeConfiguration();
    if (meshConfig->doesMeshUseSpacetree(meshName)){
      std::string spacetreeName = meshConfig->getSpacetreeName(meshName);
      meshContext->spacetree = spacetreeConfig->getSpacetree(spacetreeName);
    }
  }

  // Setup communication to server
  if (_clientMode){
    initializeClientServerCommunication();
  }
  if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode){
    initializeMasterSlaveCommunication();
  }
}

double SolverInterfaceImpl:: initialize()
{
  preciceTrace("initialize()");
  Event e("initialize", not precice::testMode);

  m2n::PointToPointCommunication::ScopedSetEventNamePrefix ssenp(
      "initialize"
      "/");

# ifndef PRECICE_NO_PETSC
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (not petscIsInitialized) {
    // Initialize Petsc if it has not already been initialized in Parallel.cpp using the precice executable.
    // This makes it possible to pass command line arguments that are consumed by Petsc when using the executable.
    PetscInitializeNoArguments();
  }
# endif

  if (_clientMode){
    preciceDebug("Request perform initializations");
    _requestManager->requestInitialize();
  }
  else {
    // Setup communication
    if (not _geometryMode){
      typedef std::map<std::string,M2NWrap>::value_type M2NPair;
      preciceInfo("initialize()", "Setting up master communication to coupling partner/s " );
      for (M2NPair& m2nPair : _m2ns) {
        m2n::M2N::SharedPointer& m2n = m2nPair.second.m2n;
        std::string localName = _accessorName;
        if (_serverMode) localName += "Server";
        std::string remoteName(m2nPair.first);
        preciceCheck(m2n.get() != NULL, "initialize()",
                     "M2N communication from " << localName << " to participant "
                     << remoteName << " could not be created! Check compile "
                     "flags used!");
        if (m2nPair.second.isRequesting){
          m2n->requestMasterConnection(remoteName, localName);
        }
        else {
          m2n->acceptMasterConnection(localName, remoteName);
        }
      }
      preciceInfo("initialize()", "Coupling partner/s are connected " );
    }

    preciceDebug("Perform initializations");
    // sort MeshContexts, provided needs to come first, and all communicated meshes must have the same order
    // on all participants
    std::sort (_accessor->usedMeshContexts().begin(), _accessor->usedMeshContexts().end(),
        []( MeshContext* lhs, const MeshContext* rhs) -> bool
        {
          if(lhs->provideMesh && not rhs->provideMesh){
            return true;
          }
          if(not lhs->provideMesh && rhs->provideMesh){
            return false;
          }
          return lhs->mesh->getName() < rhs->mesh->getName();
        } );

    for (MeshContext* meshContext : _accessor->usedMeshContexts()){
      createMeshContext(*meshContext);
    }

    if(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode){
      typedef std::map<std::string,M2NWrap>::value_type M2NPair;
      preciceInfo("initialize()", "Setting up slaves communication to coupling partner/s " );
      foreach (M2NPair& m2nPair, _m2ns){
        m2n::M2N::SharedPointer& m2n = m2nPair.second.m2n;
        std::string localName = _accessorName;
        std::string remoteName(m2nPair.first);
        preciceCheck(m2n.get() != NULL, "initialize()",
                     "Communication from " << localName << " to participant "
                     << remoteName << " could not be created! Check compile "
                     "flags used!");
        if (m2nPair.second.isRequesting){
          m2n->requestSlavesConnection(remoteName, localName);
        }
        else {
          m2n->acceptSlavesConnection(localName, remoteName);
        }
      }
    }

    std::set<action::Action::Timing> timings;
    double dt = 0.0;

    for (PtrWatchPoint& watchPoint : _accessor->watchPoints()){
      watchPoint->initialize();
    }

    // Initialize coupling state
    double time = 0.0;
    int timestep = 1;

    if (_restartMode){
      preciceInfo("initialize()", "Reading simulation state for restart");
      io::SimulationStateIO stateIO(_checkpointFileName + "_simstate.txt");
      stateIO.readState(time, timestep, _numberAdvanceCalls);
    }

    _couplingScheme->initialize(time, timestep);

    if (_restartMode){
      preciceInfo("initialize()", "Reading coupling scheme state for restart");
      //io::TXTReader txtReader(_checkpointFileName + "_cplscheme.txt");
      _couplingScheme->importState(_checkpointFileName);
    }

    dt = _couplingScheme->getNextTimestepMaxLength();

    timings.insert(action::Action::ALWAYS_POST);

    if (_couplingScheme->hasDataBeenExchanged()){
      timings.insert(action::Action::ON_EXCHANGE_POST);
      mapReadData();
    }

    performDataActions(timings, 0.0, 0.0, 0.0, dt);

    if(not utils::MasterSlave::_slaveMode){ //TODO not yet supported
      preciceDebug("Plot output...");
      for (const io::ExportContext& context : _accessor->exportContexts()){
        if (context.timestepInterval != -1){
          std::ostringstream suffix;
          suffix << _accessorName << ".init";
          exportMesh(suffix.str());
          if (context.triggerSolverPlot){
            _couplingScheme->requireAction(constants::actionPlotOutput());
          }
        }
      }
    }
    preciceInfo("initialize()", _couplingScheme->printCouplingState());
  }
  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl:: initializeData ()
{
  preciceTrace("initializeData()" );
  Event e("initializeData", not precice::testMode);

  m2n::PointToPointCommunication::ScopedSetEventNamePrefix ssenp(
      "initializeData"
      "/");

  preciceCheck(_couplingScheme->isInitialized(), "initializeData()",
               "initialize() has to be called before initializeData()");
  if (not _geometryMode){
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
    }
  }
}

double SolverInterfaceImpl:: advance
(
  double computedTimestepLength )
{
  preciceTrace1("advance()", computedTimestepLength);
  Event e("advance", not precice::testMode);

  m2n::PointToPointCommunication::ScopedSetEventNamePrefix ssenp(
      "advance"
      "/");

  preciceCheck(_couplingScheme->isInitialized(), "advance()",
               "initialize() has to be called before advance()");
  _numberAdvanceCalls++;
  preciceInfo("advance()", "Iteration #" << _numberAdvanceCalls);
  if (_clientMode){
    _requestManager->requestAdvance(computedTimestepLength);
  }
  else {

    if(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode){
      syncTimestep(computedTimestepLength);
    }

    double timestepLength = 0.0; // Length of (full) current dt
    double timestepPart = 0.0;   // Length of computed part of (full) curr. dt
    double time = 0.0;


    // Update the coupling scheme time state. Necessary to get correct remainder.
    _couplingScheme->addComputedTime(computedTimestepLength);

    if (_geometryMode){
      timestepLength = computedTimestepLength;
      timestepPart = computedTimestepLength;
    }
    else {
      //double timestepLength = 0.0;
      if (_couplingScheme->hasTimestepLength()){
        timestepLength = _couplingScheme->getTimestepLength();
      }
      else {
        timestepLength = computedTimestepLength;
      }
      timestepPart = timestepLength - _couplingScheme->getThisTimestepRemainder();
    }
    time = _couplingScheme->getTime();


    mapWrittenData();

    std::set<action::Action::Timing> timings;

    timings.insert(action::Action::ALWAYS_PRIOR);
    if (_couplingScheme->willDataBeExchanged(0.0)){
      timings.insert(action::Action::ON_EXCHANGE_PRIOR);
    }
    performDataActions(timings, time, computedTimestepLength, timestepPart, timestepLength);

    preciceDebug("Advancing coupling scheme");
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

    mapReadData();

    preciceInfo("advance()", _couplingScheme->printCouplingState());

    handleExports();

    resetWrittenData();

  }
  return _couplingScheme->getNextTimestepMaxLength();
}

void SolverInterfaceImpl:: finalize()
{
  preciceTrace("finalize()");
  preciceCheck(_couplingScheme->isInitialized(), "finalize()",
               "initialize() has to be called before finalize()");
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
          _couplingScheme->requireAction ( constants::actionPlotOutput() );
        }
      }
    }
    // Apply some final ping-pong to synch solver that run e.g. with a uni-directional coupling only
    // afterwards close connections
    typedef std::map<std::string,M2NWrap>::iterator PairIter;
    std::string ping = "ping";
    std::string pong = "pong";
    foriter ( PairIter, iter, _m2ns ){
      if( not utils::MasterSlave::_slaveMode){
        if(iter->second.isRequesting){
          iter->second.m2n->getMasterCommunication()->startSendPackage(0);
          iter->second.m2n->getMasterCommunication()->send(ping,0);
          iter->second.m2n->getMasterCommunication()->finishSendPackage();
          std::string receive = "init";
          iter->second.m2n->getMasterCommunication()->startReceivePackage(0);
          iter->second.m2n->getMasterCommunication()->receive(receive,0);
          iter->second.m2n->getMasterCommunication()->finishReceivePackage();
          assertion(receive==pong);
        }
        else{
          std::string receive = "init";
          iter->second.m2n->getMasterCommunication()->startReceivePackage(0);
          iter->second.m2n->getMasterCommunication()->receive(receive,0);
          iter->second.m2n->getMasterCommunication()->finishReceivePackage();
          assertion(receive==ping);
          iter->second.m2n->getMasterCommunication()->startSendPackage(0);
          iter->second.m2n->getMasterCommunication()->send(pong,0);
          iter->second.m2n->getMasterCommunication()->finishSendPackage();
        }
      }
      iter->second.m2n->closeConnection();
    }
  }
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    _accessor->getMasterSlaveCommunication()->closeConnection();
  }

  utils::Parallel::finalize();
}

int SolverInterfaceImpl:: getDimensions() const
{
  preciceTrace1 ( "getDimensions()", _dimensions );
  return _dimensions;
}

bool SolverInterfaceImpl:: isCouplingOngoing()
{
  preciceTrace ( "isCouplingOngoing()" );
  return _couplingScheme->isCouplingOngoing();
}

bool SolverInterfaceImpl:: isReadDataAvailable()
{
  preciceTrace ( "isReadDataAvailable()" );
  return _couplingScheme->hasDataBeenExchanged();
}

bool SolverInterfaceImpl:: isWriteDataRequired
(
  double computedTimestepLength )
{
  preciceTrace1 ( "isWriteDataRequired()", computedTimestepLength );
  return _couplingScheme->willDataBeExchanged(computedTimestepLength);
}

bool SolverInterfaceImpl:: isTimestepComplete()
{
  preciceTrace ( "isCouplingTimestepComplete()" );
  return _couplingScheme->isCouplingTimestepComplete();
}

bool SolverInterfaceImpl:: isActionRequired
(
  const std::string& action )
{
  preciceTrace2("isActionRequired()", action, _couplingScheme->isActionRequired(action));
  return _couplingScheme->isActionRequired(action);
}

void SolverInterfaceImpl:: fulfilledAction
(
  const std::string& action )
{
  preciceTrace1 ( "fulfilledAction()", action );
  if ( _clientMode ) {
    _requestManager->requestFulfilledAction(action);
  }
  _couplingScheme->performedAction(action);
}

bool SolverInterfaceImpl:: hasMesh
(
  const std::string& meshName ) const
{
  preciceTrace1 ( "hasMesh()", meshName );
  return utils::contained ( meshName, _meshIDs );
}

int SolverInterfaceImpl:: getMeshID
(
  const std::string& meshName )
{
  preciceTrace1 ( "getMeshID()", meshName );
  preciceCheck( utils::contained(meshName, _meshIDs), "getMeshID()",
                "Mesh with name \""<< meshName << "\" is not defined!" );
  return _meshIDs[meshName];
}

std::set<int> SolverInterfaceImpl:: getMeshIDs()
{
  preciceTrace ( "getMeshIDs()" );
  std::set<int> ids;
  for (const impl::MeshContext* context : _accessor->usedMeshContexts()) {
    ids.insert ( context->mesh->getID() );
  }
  return ids;
}

bool SolverInterfaceImpl:: hasData
(
  const std::string& dataName, int meshID )
{
  preciceTrace2 ( "hasData()", dataName, meshID );
  preciceCheck ( _dataIDs.find(meshID)!=_dataIDs.end(), "hasData()",
                   "No mesh with meshID \"" << meshID << "\" is defined");
  std::map<std::string,int>& sub_dataIDs =  _dataIDs[meshID];
  return sub_dataIDs.find(dataName)!= sub_dataIDs.end();
}

int SolverInterfaceImpl:: getDataID
(
  const std::string& dataName, int meshID )
{
  preciceTrace2 ( "getDataID()", dataName, meshID );
  preciceCheck ( hasData(dataName, meshID), "getDataID()",
                 "Data with name \"" << dataName << "\" is not defined on mesh with ID \""
                 << meshID << "\".");
  return _dataIDs[meshID][dataName];
}

int SolverInterfaceImpl:: inquirePosition
(
  const double*        point,
  const std::set<int>& meshIDs )
{
  preciceTrace2 ( "inquirePosition()", point, meshIDs.size() );
  using namespace precice::constants;
  int pos = positionOutsideOfGeometry();
  utils::DynVector searchPoint(_dimensions);
  for (int dim=0; dim<_dimensions; dim++) searchPoint[dim] = point[dim];
  if (_clientMode){
    pos = _requestManager->requestInquirePosition(searchPoint, meshIDs);
  }
  else {
    typedef spacetree::Spacetree Spacetree;
    std::vector<int> markedContexts(_accessor->usedMeshContexts().size());
    selectInquiryMeshIDs(meshIDs, markedContexts);
    for (int i=0; i < (int)markedContexts.size(); i++){
      MeshContext* meshContext = _accessor->usedMeshContexts()[i];
      if (markedContexts[i] == markedSkip()){
        preciceDebug("Skipping mesh " << meshContext->mesh->getName());
        continue;
      }
      int tempPos = -1;
      if (markedContexts[i] == markedQuerySpacetree()){
        assertion(meshContext->spacetree.use_count() > 0);
        tempPos = meshContext->spacetree->searchPosition(searchPoint);
      }
      else {
        assertion1(markedContexts[i] == markedQueryDirectly(), markedContexts[i]);
        query::FindClosest findClosest(searchPoint);
        findClosest(*(meshContext->mesh));
        assertion(findClosest.hasFound());
        tempPos = positionOnGeometry();
        if (tarch::la::greater(findClosest.getClosest().distance, 0.0)){
          tempPos = positionOutsideOfGeometry();
        }
        else if (tarch::la::greater(0.0, findClosest.getClosest().distance)){
          tempPos = positionInsideOfGeometry();
        }
      }

      // Union logic for multiple geometries:
      if (pos != positionInsideOfGeometry()){
        if (tempPos == positionOutsideOfGeometry()){
          if (pos != positionOnGeometry()){
            pos = tempPos; // set outside of geometry
          }
        }
        else {
          pos = tempPos; // set inside or on geometry
        }
      }
    }
  }
  preciceDebug("Return position = " << pos);
  return pos;
}

ClosestMesh SolverInterfaceImpl:: inquireClosestMesh
(
  const double*        point,
  const std::set<int>& meshIDs )
{
  preciceTrace1("inquireClosestMesh()", point);
  ClosestMesh closestMesh(_dimensions);
  utils::DynVector searchPoint(_dimensions);
  for (int dim=0; dim < _dimensions; dim++){
    searchPoint[dim] = point[dim];
  }
  if (_clientMode){
    _requestManager->requestInquireClosestMesh(searchPoint, meshIDs, closestMesh);
  }
  else {
    using namespace precice::constants;
    std::vector<int> markedContexts(_accessor->usedMeshContexts().size());
    selectInquiryMeshIDs(meshIDs, markedContexts);
    closestMesh.setPosition(positionOutsideOfGeometry());
    //foreach (MeshContext& meshContext, _accessor->usedMeshContexts()){
    for (int i=0; i < (int)markedContexts.size(); i++){
      MeshContext* meshContext = _accessor->usedMeshContexts()[i];
      if (markedContexts[i] == markedSkip()){
        preciceDebug("Skipping mesh " << meshContext->mesh->getName());
        continue;
      }
      query::FindClosest findClosest(searchPoint);
      if (markedContexts[i] == markedQuerySpacetree()){
        assertion(meshContext->spacetree.get() != NULL);
        meshContext->spacetree->searchDistance(findClosest);
      }
      else {
        assertion1(markedContexts[i] == markedQueryDirectly(), markedContexts[i]);
        findClosest(*(meshContext->mesh));
      }
      assertion(findClosest.hasFound());
      const query::ClosestElement& element = findClosest.getClosest();
      if ( element.distance > tarch::la::NUMERICAL_ZERO_DIFFERENCE &&
           closestMesh.position() == positionOutsideOfGeometry() )
      {
        if ( closestMesh.distance() > element.distance ) {
          closestMesh.setDistanceVector ( tarch::la::raw(element.vectorToElement) );
          closestMesh.meshIDs() = element.meshIDs;
        }
      }
      else if ( element.distance < - tarch::la::NUMERICAL_ZERO_DIFFERENCE ) {
        closestMesh.setPosition ( positionInsideOfGeometry() );
        if ( closestMesh.distance() > std::abs(element.distance) ) {
          closestMesh.setDistanceVector ( tarch::la::raw(element.vectorToElement) );
          closestMesh.meshIDs() = element.meshIDs;
        }
      }
      else if ( closestMesh.position() != positionInsideOfGeometry() ){
        closestMesh.setPosition ( positionOnGeometry() );
        closestMesh.setDistanceVector ( tarch::la::raw(element.vectorToElement) );
        closestMesh.meshIDs() = element.meshIDs;
      }
      //if ( _accessor->exportContext().plotNeighbors ){
      //  _exportVTKNeighbors.addNeighbors ( searchPoint, element );
      //}
    }
  }
  return closestMesh;
}

VoxelPosition SolverInterfaceImpl:: inquireVoxelPosition
(
  const double*        voxelCenter,
  const double*        voxelHalflengths,
  bool                 includeBoundaries,
  const std::set<int>& meshIDs )
{
  preciceTrace4("inquireVoxelPosition()", voxelCenter, voxelHalflengths,
                includeBoundaries, meshIDs.size());

  using namespace precice::constants;
  utils::DynVector center(_dimensions);
  utils::DynVector halflengths(_dimensions);
  for (int dim=0; dim < _dimensions; dim++){
    center[dim] = voxelCenter[dim];
    halflengths[dim] = voxelHalflengths[dim];
  }
  preciceDebug("center = " << center << ", h = " << halflengths);

  if (_clientMode){
    VoxelPosition pos;
    _requestManager->requestInquireVoxelPosition(center, halflengths, includeBoundaries, meshIDs, pos);
    return pos;
  }
  typedef spacetree::Spacetree Spacetree;
  query::FindVoxelContent::BoundaryInclusion boundaryInclude;
  boundaryInclude = includeBoundaries
                    ? query::FindVoxelContent::INCLUDE_BOUNDARY
                    : query::FindVoxelContent::EXCLUDE_BOUNDARY;

  //VoxelPosition voxelPosition;
  int pos = positionOutsideOfGeometry();
  //mesh::Group* content = new mesh::Group();
  std::vector<int> containedMeshIDs;
//  foreach (MeshContext& meshContext, _accessor->usedMeshContexts()){
//    bool skip = not utils::contained(meshContext.mesh->getID(), meshIDs);
//    skip &= not meshIDs.empty();
//    if (skip){
//      preciceDebug("Skipping mesh " << meshContext.mesh->getName());
//      continue;
//    }

  std::vector<int> markedContexts(_accessor->usedMeshContexts().size());
  selectInquiryMeshIDs(meshIDs, markedContexts);
    //foreach (MeshContext& meshContext, _accessor->usedMeshContexts()){
  for (int i=0; i < (int)markedContexts.size(); i++){
    MeshContext* meshContext = _accessor->usedMeshContexts()[i];
    if (markedContexts[i] == markedSkip()){
      preciceDebug("Skipping mesh " << meshContext->mesh->getName());
      continue;
    }
    preciceDebug("Query mesh \"" << meshContext->mesh->getName() << "\" with "
                 << meshContext->mesh->vertices().size() << " vertices");
    int oldPos = pos;
    query::FindVoxelContent findVoxel(center, halflengths, boundaryInclude);
    if (markedContexts[i] == markedQuerySpacetree()){
      assertion(meshContext->spacetree.get() != NULL);
      preciceDebug("Use spacetree for query");
      // Query first including voxel boundaries. This enables to directly
      // use cached information of spacetree cells, that do also include
      // objects on boundaries.
      pos = meshContext->spacetree->searchContent(findVoxel);

      // MERGING DISABLED!!!! CONTENT MIGHT CONTAIN DUPLICATED ELEMENTS
//      if (not findVoxel.content().empty()){
//        preciceDebug ( "Merging found content of size = " << findVoxel.content().size() );
//        mesh::Merge mergeContent;
//        mergeContent ( findVoxel.content() );
//        //findContent.content() = mergeContent.content();
//        preciceDebug ( "Merged size = " << mergeContent.content().size() );
//        content->add ( mergeContent.content() );
//      }
      //content->add(findVoxel.content());

//      if ( pos == Spacetree::ON_GEOMETRY ) {
//        query::FindVoxelContent findVoxel ( inquiryCenter, halfLengthVoxel,
//            query::FindVoxelContent::EXCLUDE_BOUNDARY );
//        findVoxel ( findVoxelInclude.content() );
//        if ( ! findVoxel.content().empty() ) {
//          content->add ( findVoxel.content() );
//        }
//      }
    }
    // The mesh does not have a spacetree
    else {
      preciceDebug("Query mesh directly");
      assertion1(markedContexts[i] == markedQueryDirectly(), markedContexts[i]);
      //query::FindVoxelContent findVoxel(center, halflengths, boundaryInclude);
      findVoxel(*meshContext->mesh);
      // If the voxel does have content
      if (not findVoxel.content().empty()){
        pos = positionOnGeometry();
        //content->add ( findVoxel.content() );
      }
      // If the voxel is empty and not inside for any other checked geometry
      else if (oldPos != positionInsideOfGeometry()){
        //preciceDebug("Query found no objects and oldpos isnt't inside");
        query::FindClosest findClosest(center);
        findClosest(*(meshContext->mesh));
        assertion(findClosest.hasFound());
        const query::ClosestElement& closest = findClosest.getClosest();
        pos = closest.distance > 0 ? positionOutsideOfGeometry()
                                   : positionInsideOfGeometry();
      }
    }

    // Retrieve mesh IDs of contained elements
    if (not findVoxel.content().empty()){
      int geoID = mesh::PropertyContainer::INDEX_GEOMETRY_ID;
      std::vector<int> tempIDs;
      std::set<int> uniqueIDs;
      for (mesh::Vertex& vertex : findVoxel.content().vertices()) {
        vertex.getProperties(geoID, tempIDs);
        for (int id : tempIDs) {
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      for (mesh::Edge& edge : findVoxel.content().edges()) {
        edge.getProperties(geoID, tempIDs);
        for (int id : tempIDs) {
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      for (mesh::Triangle& triangle : findVoxel.content().triangles()) {
        triangle.getProperties(geoID, tempIDs);
        for (int id : tempIDs) {
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      preciceDebug("Query found objects, ids.size = " << uniqueIDs.size());
      for (int id : uniqueIDs) {
        if (not utils::contained(id, containedMeshIDs)){
          containedMeshIDs.push_back(id);
        }
      }
    }

    if (oldPos == positionInsideOfGeometry()){
      preciceDebug("Since oldpos is inside, reset to inside");
      pos = positionInsideOfGeometry();
    }
    else if ((oldPos == positionOnGeometry())
             && (pos == positionOutsideOfGeometry()))
    {
      preciceDebug ( "Since old pos is on and pos is outside, reset to on" );
      pos = positionOnGeometry();
    }
    preciceDebug("pos = " << pos);
  }
  preciceDebug("Return voxel position = " << pos << ", ids.size = " << containedMeshIDs.size());
  return VoxelPosition(pos, containedMeshIDs);
}


int SolverInterfaceImpl:: getMeshVertexSize
(
  int meshID )
{
  preciceTrace1("getMeshVertexSize()", meshID);
  int size = 0;
  if (_clientMode){
    size = _requestManager->requestGetMeshVertexSize(meshID);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    assertion(context.mesh.get() != NULL);
    size = context.mesh->vertices().size();
  }
  preciceDebug("return " << size);
  return size;
}

void SolverInterfaceImpl:: resetMesh
(
  int meshID )
{
  preciceTrace1("resetMesh()", meshID);
  impl::MeshContext& context = _accessor->meshContext(meshID);
  bool hasMapping = context.fromMappingContext.mapping.use_count() > 0
            || context.toMappingContext.mapping.use_count() > 0;
  bool isStationary =
        context.fromMappingContext.timing == mapping::MappingConfiguration::INITIAL &&
            context.toMappingContext.timing == mapping::MappingConfiguration::INITIAL;

  preciceCheck(!isStationary, "resetMesh()", "A mesh with only initial mappings"
            << " must not be reseted");
  preciceCheck(hasMapping, "resetMesh()", "A mesh with no mappings"
              << " must not be reseted");

  preciceDebug ( "Clear mesh positions for mesh \"" << context.mesh->getName() << "\"" );
  context.mesh->clear ();

}

int SolverInterfaceImpl:: setMeshVertex
(
  int           meshID,
  const double* position )
{
  preciceTrace1 ( "setMeshVertex()", meshID );
  utils::DynVector internalPosition(_dimensions);
  for ( int dim=0; dim < _dimensions; dim++ ){
    internalPosition[dim] = position[dim];
  }
  preciceDebug("Position = " << internalPosition);
  int index = -1;
  if ( _clientMode ){
    index = _requestManager->requestSetMeshVertex ( meshID, internalPosition );
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    preciceDebug("MeshRequirement: " << context.meshRequirement);
    index = mesh->createVertex(internalPosition).getID();
    mesh->allocateDataValues();
  }
  return index;
}

void SolverInterfaceImpl:: setMeshVertices
(
  int     meshID,
  int     size,
  double* positions,
  int*    ids )
{
  preciceTrace2("setMeshVertices()", meshID, size);
  if (_clientMode){
    _requestManager->requestSetMeshVertices(meshID, size, positions, ids);
  }
  else { //couplingMode
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    utils::DynVector internalPosition(_dimensions);
    preciceDebug("Set positions");
    for (int i=0; i < size; i++){
      for (int dim=0; dim < _dimensions; dim++){
        internalPosition[dim] = positions[i*_dimensions + dim];
      }
      ids[i] = mesh->createVertex(internalPosition).getID();
    }
    mesh->allocateDataValues();
  }
}

void SolverInterfaceImpl:: getMeshVertices
(
  int     meshID,
  size_t  size,
  int*    ids,
  double* positions )
{
  preciceTrace2("getMeshVertices()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetMeshVertices(meshID, size, ids, positions);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    utils::DynVector internalPosition(_dimensions);
    preciceDebug("Get positions");
    assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
    for (size_t i=0; i < size; i++){
      size_t id = ids[i];
      assertion2(id < mesh->vertices().size(), mesh->vertices().size(), id);
      internalPosition = mesh->vertices()[id].getCoords();
      for (int dim=0; dim < _dimensions; dim++){
        positions[id*_dimensions + dim] = internalPosition[dim];
      }
    }
  }
}

void SolverInterfaceImpl:: getMeshVertexIDsFromPositions (
  int     meshID,
  size_t  size,
  double* positions,
  int*    ids )
{
  preciceTrace2("getMeshVertexIDsFromPositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetMeshVertexIDsFromPositions(meshID, size, positions, ids);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    mesh::PtrMesh mesh(context.mesh);
    preciceDebug("Get ids");
    utils::DynVector internalPosition(_dimensions);
    utils::DynVector position(_dimensions);
    assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
    for (size_t i=0; i < size; i++){
      for (int dim=0; dim < _dimensions; dim++){
        position[dim] = positions[i*_dimensions+dim];
      }
      size_t j=0;
      for (j=0; j < mesh->vertices().size(); j++){
        internalPosition = mesh->vertices()[j].getCoords();
        if (equals(internalPosition, position)){
          ids[i] = j;
          break;
        }
      }
      preciceCheck(j < mesh->vertices().size(), "getMeshVertexIDsFromPositions()",
                   "Position " << i << "=" << position << " unknown!");
    }
  }
}



int SolverInterfaceImpl:: setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  preciceTrace3 ( "setMeshEdge()", meshID, firstVertexID, secondVertexID );
  if (_restartMode){
    preciceDebug("Ignoring edge, since restart mode is active");
    return -1;
  }
  if ( _clientMode ){
    return _requestManager->requestSetMeshEdge ( meshID, firstVertexID, secondVertexID );
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::FULL ){
      preciceDebug("Full mesh required.");
      mesh::PtrMesh& mesh = context.mesh;
      assertion1(firstVertexID >= 0, firstVertexID);
      assertion1(secondVertexID >= 0, secondVertexID);
      assertion2(firstVertexID < (int)mesh->vertices().size(),
                 firstVertexID, mesh->vertices().size());
      assertion2(secondVertexID < (int)mesh->vertices().size(),
                 secondVertexID, mesh->vertices().size());
      mesh::Vertex& v0 = mesh->vertices()[firstVertexID];
      mesh::Vertex& v1 = mesh->vertices()[secondVertexID];
      return mesh->createEdge(v0, v1).getID ();
    }
  }
  return -1;
}

void SolverInterfaceImpl:: setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{
  preciceTrace4 ( "setMeshTriangle()", meshID, firstEdgeID,
                  secondEdgeID, thirdEdgeID );
  if (_restartMode){
    preciceDebug("Ignoring triangle, since restart mode is active");
    return;
  }
  if ( _clientMode ){
    _requestManager->requestSetMeshTriangle ( meshID, firstEdgeID, secondEdgeID, thirdEdgeID );
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.meshRequirement == mapping::Mapping::FULL ){
      mesh::PtrMesh& mesh = context.mesh;
      assertion ( firstEdgeID >= 0 );
      assertion ( secondEdgeID >= 0 );
      assertion ( thirdEdgeID >= 0 );
      assertion ( (int)mesh->edges().size() > firstEdgeID );
      assertion ( (int)mesh->edges().size() > secondEdgeID );
      assertion ( (int)mesh->edges().size() > thirdEdgeID );
      mesh::Edge& e0 = mesh->edges()[firstEdgeID];
      mesh::Edge& e1 = mesh->edges()[secondEdgeID];
      mesh::Edge& e2 = mesh->edges()[thirdEdgeID];
      mesh->createTriangle ( e0, e1, e2 );
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
  preciceTrace4("setMeshTriangleWithEdges()", meshID, firstVertexID,
                secondVertexID, thirdVertexID);
  if (_clientMode){
    _requestManager->requestSetMeshTriangleWithEdges(meshID,
                                                     firstVertexID,
                                                     secondVertexID,
                                                     thirdVertexID);
    return;
  }
  MeshContext& context = _accessor->meshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::FULL){
    mesh::PtrMesh& mesh = context.mesh;
    assertion1(firstVertexID >= 0, firstVertexID);
    assertion1(secondVertexID >= 0, secondVertexID);
    assertion1(thirdVertexID >= 0, thirdVertexID);
    assertion2((int)mesh->vertices().size() > firstVertexID,
                mesh->vertices().size(), firstVertexID);
    assertion2((int)mesh->vertices().size() > secondVertexID,
                mesh->vertices().size(), secondVertexID);
    assertion2((int)mesh->vertices().size() > thirdVertexID,
                 mesh->vertices().size(), thirdVertexID);
    mesh::Vertex* vertices[3];
    vertices[0] = &mesh->vertices()[firstVertexID];
    vertices[1] = &mesh->vertices()[secondVertexID];
    vertices[2] = &mesh->vertices()[thirdVertexID];
    mesh::Edge* edges[3];
    edges[0] = NULL;
    edges[1] = NULL;
    edges[2] = NULL;
    for (mesh::Edge& edge : mesh->edges()) {
      // Check edge 0
      bool foundEdge = edge.vertex(0).getID() == vertices[0]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[1]->getID();
      if (foundEdge){
        edges[0] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[1]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[0]->getID();
      if (foundEdge){
        edges[0] = &edge;
        continue;
      }

      // Check edge 1
      foundEdge = edge.vertex(0).getID() == vertices[1]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[2]->getID();
      if (foundEdge){
        edges[1] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[2]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[1]->getID();
      if (foundEdge){
        edges[1] = &edge;
        continue;
      }

      // Check edge 2
      foundEdge = edge.vertex(0).getID() == vertices[2]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[0]->getID();
      if (foundEdge){
        edges[2] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[0]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[2]->getID();
      if (foundEdge){
        edges[2] = &edge;
        continue;
      }
    }
    // Create missing edges
    if (edges[0] == NULL){
      edges[0] = & mesh->createEdge(*vertices[0], *vertices[1]);
    }
    if (edges[1] == NULL){
      edges[1] = & mesh->createEdge(*vertices[1], *vertices[2]);
    }
    if (edges[2] == NULL){
      edges[2] = & mesh->createEdge(*vertices[2], *vertices[0]);
    }

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
  preciceTrace5("setMeshQuad()", meshID, firstEdgeID, secondEdgeID, thirdEdgeID,
                fourthEdgeID);
  if (_restartMode){
    preciceDebug("Ignoring quad, since restart mode is active");
    return;
  }
  if (_clientMode){
    _requestManager->requestSetMeshQuad(meshID, firstEdgeID, secondEdgeID,
                                        thirdEdgeID, fourthEdgeID);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.meshRequirement == mapping::Mapping::FULL){
      mesh::PtrMesh& mesh = context.mesh;
      assertion(firstEdgeID >= 0);
      assertion(secondEdgeID >= 0);
      assertion(thirdEdgeID >= 0);
      assertion(fourthEdgeID >= 0);
      assertion((int)mesh->edges().size() > firstEdgeID);
      assertion((int)mesh->edges().size() > secondEdgeID);
      assertion((int)mesh->edges().size() > thirdEdgeID);
      assertion((int)mesh->quads().size() > fourthEdgeID);
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
  preciceTrace5("setMeshQuadWithEdges()", meshID, firstVertexID,
                secondVertexID, thirdVertexID, fourthVertexID);
  if (_clientMode){
    _requestManager->requestSetMeshQuadWithEdges(
        meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
    return;
  }
  MeshContext& context = _accessor->meshContext(meshID);
  if (context.meshRequirement == mapping::Mapping::FULL){
    mesh::PtrMesh& mesh = context.mesh;
    assertion1(firstVertexID >= 0, firstVertexID);
    assertion1(secondVertexID >= 0, secondVertexID);
    assertion1(thirdVertexID >= 0, thirdVertexID);
    assertion1(fourthVertexID >= 0, fourthVertexID);
    assertion2((int)mesh->vertices().size() > firstVertexID,
                 mesh->vertices().size(), firstVertexID);
    assertion2((int)mesh->vertices().size() > secondVertexID,
                 mesh->vertices().size(), secondVertexID);
    assertion2((int)mesh->vertices().size() > thirdVertexID,
                 mesh->vertices().size(), thirdVertexID);
    assertion2((int)mesh->vertices().size() > fourthVertexID,
                 mesh->vertices().size(), fourthVertexID);
    mesh::Vertex* vertices[4];
    vertices[0] = &mesh->vertices()[firstVertexID];
    vertices[1] = &mesh->vertices()[secondVertexID];
    vertices[2] = &mesh->vertices()[thirdVertexID];
    vertices[3] = &mesh->vertices()[fourthVertexID];
    mesh::Edge* edges[4];
    edges[0] = NULL;
    edges[1] = NULL;
    edges[2] = NULL;
    edges[3] = NULL;
    for (mesh::Edge& edge : mesh->edges()) {
      // Check edge 0
      bool foundEdge = edge.vertex(0).getID() == vertices[0]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[1]->getID();
      if ( foundEdge ){
        edges[0] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[1]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[0]->getID();
      if (foundEdge){
        edges[0] = &edge;
        continue;
      }

      // Check edge 1
      foundEdge = edge.vertex(0).getID() == vertices[1]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[2]->getID();
      if ( foundEdge ){
        edges[1] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[2]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[1]->getID();
      if ( foundEdge ){
        edges[1] = &edge;
        continue;
      }

      // Check edge 2
      foundEdge = edge.vertex(0).getID() == vertices[2]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[3]->getID();
      if ( foundEdge ){
        edges[2] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[3]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[2]->getID();
      if ( foundEdge ){
        edges[2] = &edge;
        continue;
      }

      // Check edge 3
      foundEdge = edge.vertex(0).getID() == vertices[3]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[0]->getID();
      if ( foundEdge ){
        edges[3] = &edge;
        continue;
      }
      foundEdge = edge.vertex(0).getID() == vertices[0]->getID();
      foundEdge &= edge.vertex(1).getID() == vertices[3]->getID();
      if ( foundEdge ){
        edges[3] = &edge;
        continue;
      }
    }
    // Create missing edges
    if (edges[0] == NULL){
      edges[0] = & mesh->createEdge(*vertices[0], *vertices[1]);
    }
    if (edges[1] == NULL){
      edges[1] = & mesh->createEdge(*vertices[1], *vertices[2]);
    }
    if (edges[2] == NULL){
      edges[2] = & mesh->createEdge(*vertices[2], *vertices[3]);
    }
    if (edges[3] == NULL){
      edges[3] = & mesh->createEdge(*vertices[3], *vertices[0]);
    }

    mesh->createQuad(*edges[0], *edges[1], *edges[2], *edges[3]);
  }
}

void SolverInterfaceImpl:: mapWriteDataFrom
(
  int fromMeshID )
{
  preciceTrace1("mapWriteDataFrom(int)", fromMeshID);
  if (_clientMode){
    _requestManager->requestMapWriteDataFrom(fromMeshID);
    return;
  }
  impl::MeshContext& context = _accessor->meshContext(fromMeshID);
  impl::MappingContext& mappingContext = context.fromMappingContext;
  if (mappingContext.mapping.use_count() == 0){
    preciceError("mapWriteDataFrom()", "From mesh \"" << context.mesh->getName()
                   << "\", there is no mapping defined");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    preciceDebug("Compute mapping from mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  for (impl::DataContext& context : _accessor->writeDataContexts()) {
    if (context.mesh->getID() == fromMeshID){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      assign(context.toData->values()) = 0.0;
      preciceDebug("Map data \"" << context.fromData->getName()
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
  preciceTrace1 ("mapReadDataTo(int)", toMeshID);
  if (_clientMode){
    _requestManager->requestMapReadDataTo(toMeshID);
    return;
  }
  impl::MeshContext& context = _accessor->meshContext(toMeshID);
  impl::MappingContext& mappingContext = context.toMappingContext;
  if (mappingContext.mapping.use_count() == 0){
    preciceError("mapReadDataFrom()", "From mesh \"" << context.mesh->getName()
                   << "\", there is no mapping defined!");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    preciceDebug("Compute mapping from mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  for (impl::DataContext& context : _accessor->readDataContexts()) {
    if (context.mesh->getID() == toMeshID){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      assign(context.toData->values()) = 0.0;
      preciceDebug("Map data \"" << context.fromData->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      assertion(mappingContext.mapping==context.mappingContext.mapping);
      mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.toData->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.toData->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str());
#     endif
    }
  }
  mappingContext.hasMappedData = true;
}

void SolverInterfaceImpl:: writeBlockVectorData
(
  int     fromDataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("writeBlockVectorData()", fromDataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestWriteBlockVectorData(fromDataID, size, valueIndices, values);
  }
//  else if(_slaveMode){
//      com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
//      com->send(size, 0);
//      com->send(valueIndices, size, 0);
//      com->send(values, size*_dimensions, 0);
//  }
//  else if (_masterMode){
//    preciceCheck(_accessor->isDataUsed(fromDataID), "writeBlockVectorData()",
//                 "You try to write to data that is not defined for " << _accessor->getName());
//    DataContext& context = _accessor->dataContext(fromDataID);
//    assertion(context.toData.get() != NULL);
//    utils::DynVector& valuesInternal = context.fromData->values();
//    for (int i=0; i < size; i++){
//      int offsetInternal = valueIndices[i]*_dimensions;
//      int offset = i*_dimensions;
//      for (int dim=0; dim < _dimensions; dim++){
//        assertion2(offset+dim < valuesInternal.size(),
//                   offset+dim, valuesInternal.size());
//        valuesInternal[offsetInternal + dim] = values[offset + dim];
//      }
//    }
//
//    com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
//    for(int rankSender = 0; rankSender < _accessorCommunicatorSize-1; rankSender++){
//      int slaveSize = -1;
//      com->receive(slaveSize, rankSender);
//      int* slaveIndices = new int[slaveSize];
//      com->receive(slaveIndices, slaveSize, rankSender);
//      double* slaveValues = new double[slaveSize*_dimensions];
//      com->receive(slaveValues, slaveSize*_dimensions, rankSender);
//      for (int i=0; i < slaveSize; i++){
//        int offsetInternal = slaveIndices[i]*_dimensions;
//        int offset = i*_dimensions;
//        for (int dim=0; dim < _dimensions; dim++){
//          assertion2(offset+dim < valuesInternal.size(),
//                     offset+dim, valuesInternal.size());
//          valuesInternal[offsetInternal + dim] = slaveValues[offset + dim];
//        }
//      }
//      delete[] slaveIndices;
//      delete[] slaveValues;
//    }
//
//  }
  else { //couplingMode
    preciceCheck(_accessor->isDataUsed(fromDataID), "writeBlockVectorData()",
                 "You try to write to data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(fromDataID);

    assertion(context.toData.get() != NULL);
    utils::DynVector& valuesInternal = context.fromData->values();
    for (int i=0; i < size; i++){
      int offsetInternal = valueIndices[i]*_dimensions;
      int offset = i*_dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        assertion2(offset+dim < valuesInternal.size(),
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
  preciceTrace2 ( "writeVectorData()", fromDataID, valueIndex );
# ifdef Debug
  if (_dimensions == 2) preciceDebug("value = " << tarch::la::wrap<2>(value));
  if (_dimensions == 3) preciceDebug("value = " << tarch::la::wrap<3>(value));
# endif
  preciceCheck ( valueIndex >= -1, "writeVectorData()", "Invalid value index ("
                 << valueIndex << ") when writing vector data!" );
  if (_clientMode){
    utils::DynVector valueCopy(_dimensions);
    for (int dim=0; dim < _dimensions; dim++){
      valueCopy[dim] = value[dim];
    }
    _requestManager->requestWriteVectorData(fromDataID, valueIndex, tarch::la::raw(valueCopy));
  }
  else {
    preciceCheck(_accessor->isDataUsed(fromDataID), "writeVectorData()",
             "You try to write to data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(fromDataID);
    assertion(context.toData.get() != NULL);
    utils::DynVector& values = context.fromData->values();
    assertion1(valueIndex >= 0, valueIndex);
    int offset = valueIndex * _dimensions;
    for (int dim=0; dim < _dimensions; dim++){
      values[offset+dim] = value[dim];
    }

  }
}

void SolverInterfaceImpl:: writeBlockScalarData
(
  int     fromDataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("writeBlockScalarData()", fromDataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestWriteBlockScalarData(fromDataID, size, valueIndices, values);
  }
  else {
    preciceCheck(_accessor->isDataUsed(fromDataID), "writeBlockScalarData()",
                 "You try to write to data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(fromDataID);
    assertion(context.toData.get() != NULL);
    utils::DynVector& valuesInternal = context.fromData->values();
    for (int i=0; i < size; i++){
      assertion2(i < valuesInternal.size(), i, valuesInternal.size());
      valuesInternal[valueIndices[i]] = values[i];
    }
  }
}

void SolverInterfaceImpl:: writeScalarData
(
  int    fromDataID,
  int    valueIndex,
  double value )
{
  preciceTrace3("writeScalarData()", fromDataID, valueIndex, value );
  preciceCheck(valueIndex >= -1, "writeScalarData()", "Invalid value index ("
               << valueIndex << ") when writing scalar data!");
  if (_clientMode){
    _requestManager->requestWriteScalarData(fromDataID, valueIndex, value);
  }
  else {
    preciceCheck(_accessor->isDataUsed(fromDataID), "writeScalarData()",
                 "You try to write to data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(fromDataID);
    assertion(context.toData.use_count() > 0);
    utils::DynVector& values = context.fromData->values();
    assertion1(valueIndex >= 0, valueIndex);
    values[valueIndex] = value;

  }
}

void SolverInterfaceImpl:: readBlockVectorData
(
  int     toDataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("readBlockVectorData()", toDataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestReadBlockVectorData(toDataID, size, valueIndices, values);
  }
//  else if(_slaveMode){
//    com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
//    com->send(size, 0);
//    com->send(valueIndices, size, 0);
//    com->receive(values, size*_dimensions, 0);
//  }
//  else if(_masterMode){
//    preciceCheck(_accessor->isDataUsed(toDataID), "readBlockVectorData()",
//                     "You try to read from data that is not defined for " << _accessor->getName());
//    DataContext& context = _accessor->dataContext(toDataID);
//    assertion(context.fromData.get() != NULL);
//    utils::DynVector& valuesInternal = context.toData->values();
//    for (int i=0; i < size; i++){
//      int offsetInternal = valueIndices[i] * _dimensions;
//      int offset = i * _dimensions;
//      for (int dim=0; dim < _dimensions; dim++){
//        assertion2(offsetInternal+dim < valuesInternal.size(),
//                   offsetInternal+dim, valuesInternal.size());
//        values[offset + dim] = valuesInternal[offsetInternal + dim];
//      }
//    }
//
//    com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
//    for(int rankSender = 0; rankSender < _accessorCommunicatorSize-1; rankSender++){
//      int slaveSize = -1;
//      com->receive(slaveSize, rankSender);
//      int* slaveIndices = new int[slaveSize];
//      com->receive(slaveIndices, slaveSize, rankSender);
//      double* slaveValues = new double[slaveSize*_dimensions];
//      for (int i=0; i < slaveSize; i++){
//        int offsetInternal = slaveIndices[i] * _dimensions;
//        int offset = i * _dimensions;
//        for (int dim=0; dim < _dimensions; dim++){
//          assertion2(offsetInternal+dim < valuesInternal.size(),
//                     offsetInternal+dim, valuesInternal.size());
//          slaveValues[offset + dim] = valuesInternal[offsetInternal + dim];
//        }
//      }
//      com->send(slaveValues, slaveSize*_dimensions, rankSender);
//      delete[] slaveIndices;
//      delete[] slaveValues;
//    }
//  }
  else { //couplingMode
    preciceCheck(_accessor->isDataUsed(toDataID), "readBlockVectorData()",
                 "You try to read from data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(toDataID);
    assertion(context.fromData.get() != NULL);
    utils::DynVector& valuesInternal = context.toData->values();
    for (int i=0; i < size; i++){
      int offsetInternal = valueIndices[i] * _dimensions;
      int offset = i * _dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        assertion2(offsetInternal+dim < valuesInternal.size(),
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
  double* value )
{
  preciceTrace2("readVectorData()", toDataID, valueIndex);
  preciceCheck(valueIndex >= -1, "readData(vector)", "Invalid value index ( "
               << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadVectorData(toDataID, valueIndex, value);
  }
  else {
    preciceCheck(_accessor->isDataUsed(toDataID), "readVectorData()",
                     "You try to read from data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(toDataID);
    assertion(context.fromData.use_count() > 0);
    utils::DynVector& values = context.toData->values();
    assertion1 (valueIndex >= 0, valueIndex);
    int offset = valueIndex * _dimensions;
    for (int dim=0; dim < _dimensions; dim++){
      value[dim] = values[offset + dim];
    }

  }
# ifdef Debug
  if (_dimensions == 2) preciceDebug("read value = " << tarch::la::wrap<2>(value));
  if (_dimensions == 3) preciceDebug("read value = " << tarch::la::wrap<3>(value));
# endif
}

void SolverInterfaceImpl:: readBlockScalarData
(
  int     toDataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("readBlockScalarData()", toDataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestReadBlockScalarData(toDataID, size, valueIndices, values);
  }
  else {
    preciceCheck(_accessor->isDataUsed(toDataID), "readBlockScalarData()",
                     "You try to read from data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(toDataID);
    assertion(context.fromData.get() != NULL);
    utils::DynVector& valuesInternal = context.toData->values();
    for (int i=0; i < size; i++){
      assertion2(valueIndices[i] < valuesInternal.size(),
               valueIndices[i], valuesInternal.size());
      values[i] = valuesInternal[valueIndices[i]];
    }
  }
}

void SolverInterfaceImpl:: readScalarData
(
  int     toDataID,
  int     valueIndex,
  double& value )
{
  preciceTrace3("readScalarData()", toDataID, valueIndex, value);
  preciceCheck(valueIndex >= -1, "readData(vector)", "Invalid value index ( "
               << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadScalarData(toDataID, valueIndex, value);
  }
  else {
    preciceCheck(_accessor->isDataUsed(toDataID), "readScalarData()",
                     "You try to read from data that is not defined for " << _accessor->getName());
    DataContext& context = _accessor->dataContext(toDataID);
    assertion(context.fromData.use_count() > 0);
    utils::DynVector& values = context.toData->values();
    value = values[valueIndex];

  }
  preciceDebug("Read value = " << value);
}

void SolverInterfaceImpl:: exportMesh
(
  const std::string& filenameSuffix,
  int                exportType )
{
  preciceTrace2 ( "exportMesh()", filenameSuffix, exportType );
  if ( _clientMode ){
    _requestManager->requestExportMesh ( filenameSuffix, exportType );
    return;
  }
  // Export meshes
  //const ExportContext& context = _accessor->exportContext();
  for (const io::ExportContext& context : _accessor->exportContexts()) {
    preciceDebug ( "Export type = " << exportType );
    bool exportAll = exportType == constants::exportAll();
    bool exportThis = context.exporter->getType() == exportType;
    if ( exportAll || exportThis ){
      for (MeshContext* meshContext : _accessor->usedMeshContexts()) {
        std::string name = meshContext->mesh->getName() + "-" + filenameSuffix;
        std::string filename = context.location + name;
        preciceDebug ( "Exporting mesh to file \"" << filename << "\"" );
        context.exporter->doExport ( filename, *(meshContext->mesh) );
      }
    }
    // Export spacetrees
    if (context.exportSpacetree){
      for ( MeshContext* meshContext : _accessor->usedMeshContexts()) {
        std::string name = meshContext->mesh->getName() + "-" + filenameSuffix;
        std::string filename = context.location + name + ".spacetree";
        if ( meshContext->spacetree.get() != NULL ) {
          spacetree::ExportSpacetree exportSpacetree(filename);
          exportSpacetree.doExport ( *(meshContext->spacetree) );
        }
      }
    }
  }
  // Export neighbors
  //if ( context.plotNeighbors ) {
  //  _exportVTKNeighbors.exportNeighbors ( filenameSuffix + ".neighbors" );
  //}
}


MeshHandle SolverInterfaceImpl:: getMeshHandle
(
  const std::string& meshName )
{
  preciceTrace1("getMeshHandle()", meshName);
  assertion(not _clientMode);
  for (MeshContext* context : _accessor->usedMeshContexts()){
    if (context->mesh->getName() == meshName){
      return MeshHandle(context->mesh->content());
    }
  }
  preciceError("getMeshHandle()", "Participant \"" << _accessorName
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
  preciceTrace("configureM2Ns()");
  typedef m2n::M2NConfiguration::M2NTuple M2NTuple;
  for (M2NTuple m2nTuple : config->m2ns()) {
    std::string comPartner("");
    bool isRequesting = false;
    if (boost::get<1>(m2nTuple) == _accessorName){
      comPartner = boost::get<2>(m2nTuple);
      isRequesting = true;
    }
    else if (boost::get<2>(m2nTuple) == _accessorName){
      comPartner = boost::get<1>(m2nTuple);
    }
    if (not comPartner.empty()){
      for (const impl::PtrParticipant& participant : _participants) {
        if (participant->getName() == comPartner){
          if (participant->useServer()){
            comPartner += "Server";
          }
          assertion1(not utils::contained(comPartner, _m2ns), comPartner);
          assertion(boost::get<0>(m2nTuple).use_count() > 0);
          M2NWrap m2nWrap;
          m2nWrap.m2n = boost::get<0>(m2nTuple);
          m2nWrap.isRequesting = isRequesting;
          _m2ns[comPartner] = m2nWrap;
        }
      }
    }
  }
}

void SolverInterfaceImpl:: configureSolverGeometries
(
  const m2n::M2NConfiguration::SharedPointer& m2nConfig )
{
  preciceTrace ( "configureSolverGeometries()" );
  for (MeshContext* context : _accessor->usedMeshContexts()) {
    if ( context->provideMesh ) { // Accessor provides geometry
      preciceCheck ( context->receiveMeshFrom.empty(), "configureSolverGeometries()",
                     "Participant \"" << _accessorName << "\" cannot provide "
                     << "and receive mesh " << context->mesh->getName() << "!" );
      preciceCheck ( context->geometry.use_count() == 0, "configureSolverGeometries()",
                           "Participant \"" << _accessorName << "\" cannot provide "
                           << "the geometry of mesh \"" << context->mesh->getName()
                           << " in addition to a defined geometry!" );

      bool addedReceiver = false;
      geometry::CommunicatedGeometry* comGeo = NULL;
      for (PtrParticipant receiver : _participants ) {
        for (MeshContext* receiverContext : receiver->usedMeshContexts()) {
          bool doesReceive = receiverContext->receiveMeshFrom == _accessorName;
          doesReceive &= receiverContext->mesh->getName() == context->mesh->getName();
          if ( doesReceive ){
            preciceDebug ( "   ... receiver " << receiver );
            utils::DynVector offset ( _dimensions, 0.0 );
            std::string provider ( _accessorName );

            if(!addedReceiver){
              comGeo = new geometry::CommunicatedGeometry ( offset, provider, provider,_dimensions);
              context->geometry = geometry::PtrGeometry ( comGeo );
            }
            else{
              preciceDebug ( "Further receiver added.");
            }

            // meshRequirement has to be copied from "from" to provide", since
            // mapping are only defined at "provide"
            if(receiverContext->meshRequirement > context->meshRequirement){
              context->meshRequirement = receiverContext->meshRequirement;
            }

            m2n::M2N::SharedPointer m2n =
                m2nConfig->getM2N( receiver->getName(), provider );
            comGeo->addReceiver ( receiver->getName(), m2n );
            m2n->createDistributedCommunication(context->mesh);

            addedReceiver = true;
          }
        }
      }
      if(!addedReceiver){
        preciceDebug ( "No receiver found, create SolverGeometry");
        utils::DynVector offset ( _dimensions, 0.0 );
        context->geometry = geometry::PtrGeometry (
                        new geometry::SolverGeometry ( offset) );
      }

      assertion(context->geometry.use_count() > 0);

    }
    else if ( not context->receiveMeshFrom.empty()) { // Accessor receives geometry
      preciceCheck ( not context->provideMesh, "configureSolverGeometries()",
                     "Participant \"" << _accessorName << "\" cannot provide "
                     << "and receive mesh " << context->mesh->getName() << "!" );
      utils::DynVector offset ( _dimensions, 0.0 );
      std::string receiver ( _accessorName );
      std::string provider ( context->receiveMeshFrom );
      preciceDebug ( "Receiving mesh from " << provider );
      geometry::CommunicatedGeometry * comGeo =
          new geometry::CommunicatedGeometry ( offset, receiver, provider, _dimensions );
      comGeo->setSafetyFactor(context->safetyFactor);
      m2n::M2N::SharedPointer m2n = m2nConfig->getM2N ( receiver, provider );
      comGeo->addReceiver ( receiver, m2n );
      m2n->createDistributedCommunication(context->mesh);
      preciceCheck ( context->geometry.use_count() == 0, "configureSolverGeometries()",
                     "Participant \"" << _accessorName << "\" cannot receive "
                     << "the geometry of mesh \"" << context->mesh->getName()
                     << " in addition to a defined geometry!" );
      if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
        comGeo->setBoundingFromMapping(context->fromMappingContext.mapping);
        comGeo->setBoundingToMapping(context->toMappingContext.mapping);
      }
      context->geometry = geometry::PtrGeometry ( comGeo );
    }
  }
}

void SolverInterfaceImpl:: createMeshContext
(
  MeshContext& meshContext )
{
  preciceTrace1("createMeshContext()", meshContext.mesh->getName());
  assertion ( not _clientMode );
  using boost::get;
  mesh::PtrMesh mesh = meshContext.mesh;
  geometry::PtrGeometry geometry = meshContext.geometry;
  assertion(mesh.use_count() > 0);
  std::string meshName(mesh->getName());
  if (_restartMode){
    std::string fileName("precice_checkpoint_" + _accessorName + "_" + meshName);
    geometry::ImportGeometry* importGeo = new geometry::ImportGeometry (
        utils::DynVector(_dimensions, 0.0), fileName,
        geometry::ImportGeometry::VRML_1_FILE, true, not meshContext.provideMesh);
    geometry = geometry::PtrGeometry ( importGeo );
  }
  else if ( (not _geometryMode) && (geometry.use_count() > 0) ){
    utils::DynVector offset(geometry->getOffset());
    offset += meshContext.localOffset;
    preciceDebug("Adding local offset = " << meshContext.localOffset
                 << " to mesh " << mesh->getName());
    geometry->setOffset(offset);
  }

  assertion(not (_geometryMode && (geometry.use_count() == 0)));
  if (geometry.use_count() > 0){
    geometry->create(*mesh);
    preciceDebug("Created geometry \"" << meshName
                 << "\" with # vertices = " << mesh->vertices().size());
  }

  // Create spacetree for the geometry, if configured so
  if (meshContext.spacetree.use_count() > 0){
    preciceCheck(_geometryMode, "createMeshContext()",
                 "Creating spacetree in coupling mode!");
    meshContext.spacetree->addMesh(mesh);
  }
}

void SolverInterfaceImpl:: mapWrittenData()
{
  preciceTrace("mapWrittenData()");
  using namespace mapping;
  MappingConfiguration::Timing timing;
  // Compute mappings
  for (impl::MappingContext& context : _accessor->writeMappingContexts()) {
    timing = context.timing;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    bool hasComputed = context.mapping->hasComputedMapping();
    bool isNotEmpty = not _accessor->meshContext(context.fromMeshID).mesh->vertices().empty();
    if (rightTime && not hasComputed && isNotEmpty){
      preciceDebug("Compute write mapping from mesh \""
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
    bool hasMapping = context.mappingContext.mapping.get() != NULL;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    bool hasMapped = context.mappingContext.hasMappedData;
    bool isNotEmpty = context.fromData->values().size()>0;
    if (hasMapping && rightTime && (not hasMapped) && isNotEmpty){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      preciceDebug("Map data \"" << context.fromData->getName()
                   << "\" from mesh \"" << context.mesh->getName() << "\"");
      assign(context.toData->values()) = 0.0;
      preciceDebug("Map from dataID " << inDataID << " to dataID: " << outDataID);
      context.mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.fromData->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.toData->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str() );
#     endif
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
  preciceTrace("mapReadData()");
  mapping::MappingConfiguration::Timing timing;
  // Compute mappings
  for (impl::MappingContext& context : _accessor->readMappingContexts()) {
    timing = context.timing;
    bool mapNow = timing == mapping::MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == mapping::MappingConfiguration::INITIAL;
    bool hasComputed = context.mapping->hasComputedMapping();
    bool isNotEmpty = not _accessor->meshContext(context.toMeshID).mesh->vertices().empty();
    if (mapNow && not hasComputed && isNotEmpty){
      preciceDebug("Compute read mapping from mesh \""
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
    bool hasMapping = context.mappingContext.mapping.get() != NULL;
    bool hasMapped = context.mappingContext.hasMappedData;
    bool isNotEmpty = context.toData->values().size()>0;
    if (mapNow && hasMapping && (not hasMapped) && isNotEmpty){
      int inDataID = context.fromData->getID();
      int outDataID = context.toData->getID();
      assign(context.toData->values()) = 0.0;
      preciceDebug("Map read data \"" << context.fromData->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      context.mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.toData->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.toData->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str());
#     endif
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
  preciceTrace("performDataActions()");
  assertion(not _clientMode);
  for (action::PtrAction& action : _accessor->actions()) {
    if (timings.find(action->getTiming()) != timings.end()){
      action->performAction(time, dt, partFullDt, fullDt);
    }
  }
}

void SolverInterfaceImpl:: handleExports()
{
  preciceTrace("handleExports()");
  assertion(not _clientMode);
  //timesteps was already incremented before
  int timesteps = _couplingScheme->getTimesteps()-1;

  if(not utils::MasterSlave::_slaveMode){ //TODO  not yet supported
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
              _couplingScheme->requireAction(constants::actionPlotOutput());
            }
          }
        }
      }
    }
  }

  if (_couplingScheme->isCouplingTimestepComplete()){
    // Export watch point data
    for (PtrWatchPoint watchPoint : _accessor->watchPoints()) {
      watchPoint->exportPointData(_couplingScheme->getTime());
    }

    if(not utils::MasterSlave::_slaveMode){ //TODO not yet supported
      // Checkpointing
      int checkpointingInterval = _couplingScheme->getCheckpointTimestepInterval();
      preciceCheck(not  (utils::MasterSlave::_masterMode && checkpointingInterval!=-1) ,
                    "handleExports()","Checkpointing for a Master is not yet supported");
      if ((checkpointingInterval != -1) && (timesteps % checkpointingInterval == 0)){
        preciceDebug("Set require checkpoint");
        _couplingScheme->requireAction(constants::actionWriteSimulationCheckpoint());
        for (const MeshContext* meshContext : _accessor->usedMeshContexts()) {
          io::ExportVRML exportVRML(false);
          std::string filename("precice_checkpoint_" + _accessorName
                               + "_" + meshContext->mesh->getName());
          exportVRML.doExportCheckpoint(filename, *(meshContext->mesh));
        }
        io::SimulationStateIO exportState(_checkpointFileName + "_simstate.txt");

        exportState.writeState(_couplingScheme->getTime(),_couplingScheme->getTimesteps(), _numberAdvanceCalls);
        //io::TXTWriter exportCouplingSchemeState(_checkpointFileName + "_cplscheme.txt");
        _couplingScheme->exportState(_checkpointFileName);
      }
    }
  }
}

void SolverInterfaceImpl:: resetWrittenData()
{
  preciceTrace("resetWrittenData()");
  for (DataContext& context : _accessor->writeDataContexts()) {
    assign(context.fromData->values()) = 0.0;
    if (context.toData != context.fromData){
      assign(context.toData->values()) = 0.0;
    }
  }
//  if ( _accessor->exportContext().plotNeighbors ){
//    _exportVTKNeighbors.resetElements ();
//  }
}

//void SolverInterfaceImpl:: resetDataIndices()
//{
//  preciceTrace ( "resetDataIndices()" );
//  foreach ( DataContext & context, _accessor->writeDataContexts() ){
//    context.indexCursor = 0;
//  }
//  foreach ( DataContext & context, _accessor->readDataContexts() ){
//    context.indexCursor = 0;
//  }
//}

PtrParticipant SolverInterfaceImpl:: determineAccessingParticipant
(
   const config::SolverInterfaceConfiguration& config )
{
  config::PtrParticipantConfiguration partConfig =
      config.getParticipantConfiguration ();
  for (const PtrParticipant& participant : partConfig->getParticipants()) {
    if ( participant->getName() == _accessorName ) {
      return participant;
    }
  }
  preciceError ( "determineAccessingParticipant()",
                 "Accessing participant \"" << _accessorName << "\" is not defined"
                 << " in configuration!" );
}

void SolverInterfaceImpl:: selectInquiryMeshIDs
(
  const std::set<int>& meshIDs,
  std::vector<int>&    markedMeshContexts ) const
{
  preciceTrace1("selectInquiryMeshIDs()", meshIDs.size());
  assertion2(markedMeshContexts.size() == _accessor->usedMeshContexts().size(),
             markedMeshContexts.size(), _accessor->usedMeshContexts().size());

  if (meshIDs.empty()){ // All mesh IDs are used in inquiry
    for (int i=0; i < (int)markedMeshContexts.size(); i++){
      const MeshContext* context = _accessor->usedMeshContexts()[i];
      if (context->spacetree.get() == NULL){
        markedMeshContexts[i] = markedQueryDirectly();
      }
      else if (context->mesh->getID() == context->spacetree->meshes().front()->getID()){
        markedMeshContexts[i] = markedQuerySpacetree();
      }
      else {
        markedMeshContexts[i] = markedSkip();
      }
    }
  }
  else {
    for (int i=0; i < (int)markedMeshContexts.size(); i++){
      const MeshContext* context = _accessor->usedMeshContexts()[i];
      if (utils::contained(context->mesh->getID(), meshIDs)){
        if (context->spacetree.get() == NULL){
          markedMeshContexts[i] = markedQueryDirectly();
        }
        else {
          bool allSpacetreeMeshesAreInquired = true;
          for (const mesh::PtrMesh& mesh : context->spacetree->meshes()) {
            if (not utils::contained(mesh->getID(), meshIDs)){
              allSpacetreeMeshesAreInquired = false;
              break;
            }
          }
          if (allSpacetreeMeshesAreInquired){
            bool isFirst = context->mesh->getID()
                           == context->spacetree->meshes().front()->getID();
            if (isFirst){
              markedMeshContexts[i] = markedQuerySpacetree();
            }
            else {
              // Not selected, since already covered by query of first spacetree
              // mesh.
              markedMeshContexts[i] = markedSkip();
            }
          }
          else {
            markedMeshContexts[i] = markedQueryDirectly();
          }
        }
      }
      else {
        markedMeshContexts[i] = markedSkip();
      }
    }
  }
}

void SolverInterfaceImpl:: initializeClientServerCommunication()
{
  preciceTrace ( "initializeClientServerCom.()" );
  com::Communication::SharedPointer com = _accessor->getClientServerCommunication();
  assertion(com.get() != NULL);
  if ( _serverMode ){
    preciceInfo ( "initializeClientServerCom.()", "Setting up communication to client" );
    com->acceptConnection ( _accessorName + "Server", _accessorName,
                            _accessorProcessRank, _accessorCommunicatorSize );
  }
  else {
    preciceInfo ( "initializeClientServerCom.()", "Setting up communication to server" );
    com->requestConnection( _accessorName + "Server", _accessorName,
                            _accessorProcessRank, _accessorCommunicatorSize );
  }
}

void SolverInterfaceImpl:: initializeMasterSlaveCommunication()
{
  preciceTrace ( "initializeMasterSlaveCom.()" );
  com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
  assertion(com.get() != NULL);
  utils::MasterSlave::_communication = com;
  //slaves create new communicator with ranks 0 to size-2
  //therefore, the master uses a rankOffset and the slaves have to call request
  // with that offset
  int rankOffset = 1;
  if ( utils::MasterSlave::_masterMode ){
    preciceInfo ( "initializeMasterSlaveCom.()", "Setting up communication to slaves" );
    com->acceptConnection ( _accessorName + "Master", _accessorName,
                            _accessorProcessRank, 1);
    com->setRankOffset(rankOffset);
  }
  else {
    assertion(utils::MasterSlave::_slaveMode);
    com->requestConnection( _accessorName + "Master", _accessorName,
                            _accessorProcessRank-rankOffset, _accessorCommunicatorSize-rankOffset );
  }
}

void SolverInterfaceImpl:: syncTimestep(double computedTimestepLength)
{
  assertion(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode);
  com::Communication::SharedPointer com = _accessor->getMasterSlaveCommunication();
  if(utils::MasterSlave::_slaveMode){
    com->send(computedTimestepLength, 0);
  }
  else if(utils::MasterSlave::_masterMode){
    for(int rankSlave = 1; rankSlave < _accessorCommunicatorSize; rankSlave++){
      double dt;
      com->receive(dt, rankSlave);
      preciceCheck(tarch::la::equals(dt, computedTimestepLength), "advance()",
                 "Ambiguous timestep length when calling request advance from "
                 << "several processes!");
    }
  }
}


}} // namespace precice, impl

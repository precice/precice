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
#include "com/config/CommunicationConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "geometry/Geometry.hpp"
#include "geometry/ImportGeometry.hpp"
#include "geometry/CommunicatedGeometry.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "utils/String.hpp"
#include "mapping/Mapping.hpp"
#include <set>
#include <limits>
#include <cstring>
#include "boost/tuple/tuple.hpp"

namespace precice {
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
  _communications(),
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
}

SolverInterfaceImpl:: ~SolverInterfaceImpl()
{
  if (_requestManager != NULL){
    delete _requestManager;
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
  _clientMode = (not _serverMode) && _accessor->useServer();
  _participants = config.getParticipantConfiguration()->getParticipants();
  configureCommunications(config.getCommunicationConfiguration());

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
  }
  else if (not _clientMode){
    preciceInfo("configure()", "[PRECICE] Run in coupling mode");
    preciceCheck(_participants.size() > 1,
                 "configure()", "At least two participants need to be defined!");
    configureCommunicatedGeometries(config.getCommunicationConfiguration());
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
    com::PtrCommunication com = _accessor->getClientServerCommunication();
    assertion(com.get() != NULL);
    _requestManager = new RequestManager(_geometryMode, *this, com, _couplingScheme);
  }

  // Add meshIDs, data IDs, and spacetrees
  foreach (MeshContext& meshContext, _accessor->usedMeshContexts()){
    const mesh::PtrMesh& mesh = meshContext.mesh;
    std::pair<std::string,int> nameID;
    foreach (nameID, mesh->getNameIDPairs()){
      assertion(not utils::contained(nameID.first, _meshIDs));
      _meshIDs[nameID.first] = nameID.second;
    }
    foreach (const mesh::PtrData& data, mesh->data()){
      if (not utils::contained(data->getName(), _dataIDs)){
        _dataIDs[data->getName()] = data->getID();
      }
      else {
        preciceDebug("DataID: " << data->getID() << " _dataIDs: "
                     << _dataIDs[data->getName()]);
        assertion2(_dataIDs[data->getName()] == data->getID(),
                   _dataIDs[data->getName()], data->getID());
      }
    }
    std::string meshName = mesh->getName();
    mesh::PtrMeshConfiguration meshConfig = config.getMeshConfiguration();
    spacetree::PtrSpacetreeConfiguration spacetreeConfig = config.getSpacetreeConfiguration();
    if (meshConfig->doesMeshUseSpacetree(meshName)){
      std::string spacetreeName = meshConfig->getSpacetreeName(meshName);
      meshContext.spacetree = spacetreeConfig->getSpacetree(spacetreeName);
    }
  }

  // Setup communication to server
  if (_clientMode){
    initializeClientServerCommunication();
  }
}

double SolverInterfaceImpl:: initialize()
{
  preciceTrace("initialize()");

  // Perform initializations
  if (_clientMode){
    preciceDebug("Request perform initializations");
    _requestManager->requestInitialize();
  }
  else {
    // Setup communication
    if (not _geometryMode){
      preciceInfo("initialize()", "Setting up communication to coupling partner/s");
      typedef std::map<std::string,Communication>::value_type ComPair;
      foreach (ComPair& comPair, _communications){
        com::PtrCommunication& communication = comPair.second.communication;
        std::string localName = _accessorName;
        if (_serverMode) localName += "Server";
        std::string remoteName(comPair.first);
        preciceCheck(communication.get() != NULL, "initialize()",
                     "Communication from " << localName << " to participant "
                     << remoteName << " could not be created! Check compile "
                     "flags used!");
        if (comPair.second.isRequesting){
          communication->requestConnection(remoteName, localName,
              _accessorProcessRank, _accessorCommunicatorSize);
        }
        else {
          communication->acceptConnection(localName, remoteName,
              _accessorProcessRank, _accessorCommunicatorSize);
        }
      }
    }

    preciceDebug("Perform initializations");
    foreach (MeshContext& meshContext, _accessor->usedMeshContexts()){
      createMeshContext(meshContext);
    }
    foreach (PtrWatchPoint& watchPoint, _accessor->watchPoints()){
      watchPoint->initialize();
    }

    // Initialize coupling state
    double time = 0.0;
    int timestep = 0;
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
    double dt = _couplingScheme->getNextTimestepMaxLength();
    std::set<action::Action::Timing> timings;
    timings.insert(action::Action::ALWAYS_POST);
    if (_couplingScheme->hasDataBeenExchanged()){
      timings.insert(action::Action::ON_EXCHANGE_POST);
      mapReadData();
    }
    performDataActions(timings, 0.0, 0.0, 0.0, dt);
    preciceDebug("Plot output...");
    foreach (const io::ExportContext& context, _accessor->exportContexts()){
      if (context.timestepInterval != -1){
        std::ostringstream suffix;
        suffix << _accessorName << ".init";
        exportMesh(suffix.str());
        if (context.triggerSolverPlot){
          _couplingScheme->requireAction(constants::actionPlotOutput());
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
  preciceCheck(_couplingScheme->isInitialized(), "advance()",
               "initialize() has to be called before advance()");
  _numberAdvanceCalls++;
  if (_clientMode){
    _requestManager->requestAdvance(computedTimestepLength);
  }
  else {
    // Update the coupling scheme time state. Necessary to get correct remainder.
    _couplingScheme->addComputedTime(computedTimestepLength);

    double timestepLength = 0.0; // Length of (full) current dt
    double timestepPart = 0.0;   // Length of computed part of (full) curr. dt
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
    double time = _couplingScheme->getTime();
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
    foreach (const io::ExportContext& context, _accessor->exportContexts()){
      if ( context.timestepInterval != -1 ){
        std::ostringstream suffix;
        suffix << _accessorName << ".final";
        exportMesh ( suffix.str() );
        if ( context.triggerSolverPlot ) {
          _couplingScheme->requireAction ( constants::actionPlotOutput() );
        }
      }
    }
    typedef std::map<std::string,Communication>::iterator PairIter;
    foriter ( PairIter, iter, _communications ){
      iter->second.communication->closeConnection();
    }
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
  foreach ( const impl::MeshContext& context, _accessor->usedMeshContexts() ){
    ids.insert ( context.mesh->getID() );
  }
  return ids;
}

bool SolverInterfaceImpl:: hasData
(
  const std::string& dataName ) const
{
  preciceTrace1 ( "hasData()", dataName );
  return utils::contained ( dataName, _dataIDs );
}

int SolverInterfaceImpl:: getDataID
(
  const std::string& dataName )
{
  preciceTrace1 ( "getDataID()", dataName );
  preciceCheck ( utils::contained(dataName, _dataIDs), "getDataID()",
                 "Data with name \"" << dataName << "\" is not defined!" );
  return _dataIDs[dataName];
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
      MeshContext& meshContext = _accessor->usedMeshContexts()[i];
      if (markedContexts[i] == markedSkip()){
        preciceDebug("Skipping mesh " << meshContext.mesh->getName());
        continue;
      }
      int tempPos = -1;
      if (markedContexts[i] == markedQuerySpacetree()){
        assertion(meshContext.spacetree.use_count() > 0);
        tempPos = meshContext.spacetree->searchPosition(searchPoint);
      }
      else {
        assertion1(markedContexts[i] == markedQueryDirectly(), markedContexts[i]);
        query::FindClosest findClosest(searchPoint);
        findClosest(*(meshContext.mesh));
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
      MeshContext& meshContext = _accessor->usedMeshContexts()[i];
      if (markedContexts[i] == markedSkip()){
        preciceDebug("Skipping mesh " << meshContext.mesh->getName());
        continue;
      }
      query::FindClosest findClosest(searchPoint);
      if (markedContexts[i] == markedQuerySpacetree()){
        assertion(meshContext.spacetree.get() != NULL);
        meshContext.spacetree->searchDistance(findClosest);
      }
      else {
        assertion1(markedContexts[i] == markedQueryDirectly(), markedContexts[i]);
        findClosest(*(meshContext.mesh));
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
    MeshContext& meshContext = _accessor->usedMeshContexts()[i];
    if (markedContexts[i] == markedSkip()){
      preciceDebug("Skipping mesh " << meshContext.mesh->getName());
      continue;
    }
    preciceDebug("Query mesh \"" << meshContext.mesh->getName() << "\" with "
                 << meshContext.mesh->vertices().size() << " vertices");
    int oldPos = pos;
    query::FindVoxelContent findVoxel(center, halflengths, boundaryInclude);
    if (markedContexts[i] == markedQuerySpacetree()){
      assertion(meshContext.spacetree.get() != NULL);
      preciceDebug("Use spacetree for query");
      // Query first including voxel boundaries. This enables to directly
      // use cached information of spacetree cells, that do also include
      // objects on boundaries.
      pos = meshContext.spacetree->searchContent(findVoxel);

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
      findVoxel(*meshContext.mesh);
      // If the voxel does have content
      if (not findVoxel.content().empty()){
        pos = positionOnGeometry();
        //content->add ( findVoxel.content() );
      }
      // If the voxel is empty and not inside for any other checked geometry
      else if (oldPos != positionInsideOfGeometry()){
        //preciceDebug("Query found no objects and oldpos isnt't inside");
        query::FindClosest findClosest(center);
        findClosest(*(meshContext.mesh));
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
      foreach (mesh::Vertex& vertex, findVoxel.content().vertices()){
        vertex.getProperties(geoID, tempIDs);
        foreach (int id, tempIDs){
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      foreach (mesh::Edge& edge, findVoxel.content().edges()){
        edge.getProperties(geoID, tempIDs);
        foreach (int id, tempIDs){
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      foreach (mesh::Triangle& triangle, findVoxel.content().triangles()){
        triangle.getProperties(geoID, tempIDs);
        foreach (int id, tempIDs){
          uniqueIDs.insert(id);
        }
        tempIDs.clear();
      }
      preciceDebug("Query found objects, ids.size = " << uniqueIDs.size());
      foreach (int id, uniqueIDs){
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

int SolverInterfaceImpl:: setMeshVertex
(
  int           meshID,
  const double* position )
{
  preciceTrace1("setMeshVertex()", meshID);
  int index = -1;
  if (_restartMode){
    preciceDebug("Ignoring vertex, since restart mode is active");
    return index;
  }
  utils::DynVector internalPosition(_dimensions);
  for ( int dim=0; dim < _dimensions; dim++ ){
    internalPosition[dim] = position[dim];
  }
  preciceDebug( "position = " << internalPosition );
  if ( _clientMode ){
    index = _requestManager->requestSetMeshVertex ( meshID, internalPosition );
  }
  else {
    MeshContext & context = _accessor->meshContext ( meshID );
    if (context.meshRequirement == mapping::Mapping::FULL){
      preciceDebug("Set mesh vertex");
      assertion(context.writeMappingContext.mapping.use_count() == 0);
      assertion(context.readMappingContext.mapping.use_count() == 0);
      index = context.mesh->createVertex(internalPosition).getID();
      context.mesh->allocateDataValues();
    }
    preciceDebug("Vertex index = " << index);
  }
  return index;
}

void SolverInterfaceImpl:: resetWritePositions
(
  int meshID )
{
  preciceTrace1("resetWritePositions()", meshID);
  impl::MeshContext& context = _accessor->meshContext(meshID);
  bool hasMapping = context.writeMappingContext.mapping.use_count() > 0;
  bool isIncremental = context.writeMappingContext.timing
                       == mapping::MappingConfiguration::INCREMENTAL;
  bool isStationary = context.writeMappingContext.timing
                      == mapping::MappingConfiguration::INITIAL;
  if ( hasMapping && (not isIncremental) && (not isStationary) ) {
    preciceDebug ( "Clear write mesh positions for mesh \""
                   << context.mesh->getName() << "\"" );
    context.writeMappingContext.localMesh->clear ();
  }
}

void SolverInterfaceImpl:: resetReadPositions
(
  int meshID )
{
  preciceTrace1 ( "resetReadPositions()", meshID );
  impl::MeshContext& context = _accessor->meshContext(meshID);
  bool hasMapping = context.readMappingContext.mapping.use_count() > 0;
  bool isIncremental = context.readMappingContext.timing
                       == mapping::MappingConfiguration::INCREMENTAL;
  bool isStationary = context.readMappingContext.timing
                      == mapping::MappingConfiguration::INITIAL;
  if ( hasMapping && (not isIncremental) && (not isStationary)  ) {
    preciceDebug ( "Clear read mesh positions for mesh \""
                   << context.mesh->getName() << "\"" );
    context.readMappingContext.localMesh->clear ();
  }
}

int SolverInterfaceImpl:: setWritePosition
(
  int           meshID,
  const double* position )
{
  preciceTrace1 ( "setWritePosition()", meshID );
  utils::DynVector internalPosition(_dimensions);
  for ( int dim=0; dim < _dimensions; dim++ ){
    internalPosition[dim] = position[dim];
  }
  preciceDebug("Position = " << internalPosition);
  int index = -1;
  if ( _clientMode ){
    index = _requestManager->requestSetWritePosition ( meshID, internalPosition );
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.writeMappingContext.mapping.get() == NULL){
      preciceDebug("No write position required");
    }
    else {
      mesh::PtrMesh mesh(context.writeMappingContext.localMesh);
      if (context.writeMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Set temporary write position");
        assertion(mesh->vertices().size() == 1);
        mesh->vertices()[0].setCoords(internalPosition);
        context.writeMappingContext.mapping->computeMapping();
        index = 0;
      }
      else {
        preciceDebug("Set write position");
        index = mesh->createVertex(internalPosition).getID();
        mesh->allocateDataValues();
      }
    }
  }
  return index;
}

void SolverInterfaceImpl:: setWritePositions
(
  int     meshID,
  int     size,
  double* positions,
  int*    ids )
{
  preciceTrace2("setWritePositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestSetWritePositions(meshID, size, positions, ids);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.writeMappingContext.mapping.get() == NULL){
      preciceDebug("No write position required");
    }
    else {
      mesh::PtrMesh mesh(context.writeMappingContext.localMesh);
      utils::DynVector internalPosition(_dimensions);
      if (context.writeMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Set temporary write position");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        for ( int dim=0; dim < _dimensions; dim++ ){
          internalPosition[dim] = positions[dim];
        }
        mesh->vertices()[0].setCoords(internalPosition);
        context.writeMappingContext.mapping->computeMapping();
        ids[0] = 0;
      }
      else {
        preciceDebug("Set write positions");
        for (int i=0; i < size; i++){
          for (int dim=0; dim < _dimensions; dim++){
            internalPosition[dim] = positions[i*_dimensions + dim];
          }
          ids[i] = mesh->createVertex(internalPosition).getID();
        }
        mesh->allocateDataValues();
      }
    }
  }
}

void SolverInterfaceImpl:: getWritePositions
(
  int     meshID,
  int     size,
  int*    ids,
  double* positions )
{
  preciceTrace2("getWritePositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetWritePositions(meshID, size, ids, positions);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.writeMappingContext.mapping.get() == NULL){
      preciceWarning("getWritePositions()", "No write positions available!");
    }
    else {
      mesh::PtrMesh mesh(context.writeMappingContext.localMesh);
      utils::DynVector internalPosition(_dimensions);
      if (context.writeMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Get temporary write position");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        internalPosition = mesh->vertices()[0].getCoords();
        for (int dim=0; dim < _dimensions; dim++){
          positions[dim] = internalPosition[dim];
        }
      }
      else {
        preciceDebug("Get write positions");
        assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
        for (int i=0; i < size; i++){
          int id = ids[i];
          assertion2(id < mesh->vertices().size(), mesh->vertices().size(), id);
          internalPosition = mesh->vertices()[id].getCoords();
          for (int dim=0; dim < _dimensions; dim++){
            positions[id*_dimensions + dim] = internalPosition[dim];
          }
        }
      }
    }
  }
}

void SolverInterfaceImpl:: getWriteIDsFromPositions (
  int     meshID,
  int     size,
  double* positions,
  int*    ids )
{
  preciceTrace2("getWriteIDsFromPositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetWriteIDsFromPositions(meshID, size, positions, ids);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.writeMappingContext.mapping.get() == NULL){
      preciceWarning("getWriteIDsFromPositions()", "No write ids available!");
    }
    else {
      mesh::PtrMesh mesh(context.writeMappingContext.localMesh);
      if (context.writeMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Get temporary write id --> 0");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        ids[0] = 0;
      }
      else {
        preciceDebug("Get write ids");
        utils::DynVector internalPosition(_dimensions);
        utils::DynVector position(_dimensions);
        assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
        for (int i=0; i < size; i++){
          for (int dim=0; dim < _dimensions; dim++){
            position[dim] = positions[i*_dimensions+dim];
          }
          int j=0;
          for (; j < mesh->vertices().size(); j++){
            internalPosition = mesh->vertices()[j].getCoords();
            if (equals(internalPosition, position)){
              ids[i] = j;
              break;
            }
          }
          preciceCheck(j < mesh->vertices().size(), "getWriteIDsFromPositions()",
                       "Position " << i << "=" << position << " unknown!");
        }
      }
    }
  }
}


int SolverInterfaceImpl:: getWriteNodesSize
(
  int meshID )
{
  preciceTrace1("getWriteNodesSize()", meshID);
  if (_clientMode){
    return _requestManager->requestGetWriteNodesSize(meshID);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.writeMappingContext.mapping.get() == NULL){
      return 0;
    }
    else {
      mesh::PtrMesh mesh(context.writeMappingContext.localMesh);
      return mesh->vertices().size();
    }
  }
}


int SolverInterfaceImpl:: setReadPosition
(
  int           meshID,
  const double* position )
{
  preciceTrace1("setReadPosition()", meshID);
  utils::DynVector internalPosition(_dimensions);
  for ( int dim=0; dim < _dimensions; dim++ ){
    internalPosition[dim] = position[dim];
  }
  preciceDebug("Position = " << internalPosition);
  int index = -1;
  if ( _clientMode ){
    index = _requestManager->requestSetReadPosition(meshID, internalPosition);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if ( context.readMappingContext.mapping.get() == NULL){
      preciceDebug("No read position required");
      index = -1;
    }
    else {
      mesh::PtrMesh mesh(context.readMappingContext.localMesh);
      if ( context.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL ){
        preciceDebug("Set temporary read position");
        assertion1(mesh->vertices().size() == 1, mesh->vertices().size());
        mesh->vertices()[0].setCoords(internalPosition);
        context.readMappingContext.mapping->computeMapping();
        index = 0;
      }
      else {
        preciceDebug("Set read position");
        index = mesh->createVertex(internalPosition).getID();
        mesh->allocateDataValues();
      }
    }
  }
  return index;
}

void SolverInterfaceImpl:: setReadPositions
(
  int     meshID,
  int     size,
  double* positions,
  int*    ids )
{
  preciceTrace2("setReadPositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestSetReadPositions(meshID, size, positions, ids);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.readMappingContext.mapping.get() == NULL){
      preciceWarning("setReadPositions()", "No read position required!");
    }
    else {
      mesh::PtrMesh mesh(context.readMappingContext.localMesh);
      utils::DynVector internalPosition(_dimensions);
      if (context.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Set temporary read position");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        for ( int dim=0; dim < _dimensions; dim++ ){
          internalPosition[dim] = positions[dim];
        }
        mesh->vertices()[0].setCoords(internalPosition);
        context.readMappingContext.mapping->computeMapping();
        ids[0] = 0;
      }
      else {
        preciceDebug("Set read positions");
        for (int i=0; i < size; i++){
          for (int dim=0; dim < _dimensions; dim++){
            internalPosition[dim] = positions[i*_dimensions + dim];
          }
          ids[i] = mesh->createVertex(internalPosition).getID();
        }
        mesh->allocateDataValues();
      }
    }
  }
}

void SolverInterfaceImpl:: getReadPositions
(
  int     meshID,
  int     size,
  int*    ids,
  double* positions )
{
  preciceTrace2("getReadPositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetReadPositions(meshID, size, ids, positions);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.readMappingContext.mapping.get() == NULL){
      preciceWarning("getReadPositions()", "No read positions available!");
    }
    else {
      mesh::PtrMesh mesh(context.readMappingContext.localMesh);
      utils::DynVector internalPosition(_dimensions);
      if (context.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Get temporary read position");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        internalPosition = mesh->vertices()[0].getCoords();
        for (int dim=0; dim < _dimensions; dim++){
          positions[dim] = internalPosition[dim];
        }
      }
      else {
        preciceDebug("Get read positions");
        assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
        for (int i=0; i < size; i++){
          int id = ids[i];
          assertion2(id < mesh->vertices().size(), mesh->vertices().size(), id);
          internalPosition = mesh->vertices()[id].getCoords();
          for (int dim=0; dim < _dimensions; dim++){
            positions[id*_dimensions + dim] = internalPosition[dim];
          }
        }
      }
    }
  }
}

void SolverInterfaceImpl:: getReadIDsFromPositions (
  int     meshID,
  int     size,
  double* positions,
  int*    ids )
{
  preciceTrace2("getReadIDsFromPositions()", meshID, size);
  if (_clientMode){
    _requestManager->requestGetReadIDsFromPositions(meshID, size, positions, ids);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.readMappingContext.mapping.get() == NULL){
      preciceWarning("getReadIDsFromPositions()", "No read ids available!");
    }
    else {
      mesh::PtrMesh mesh(context.readMappingContext.localMesh);
      if (context.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
        preciceDebug("Get temporary read id --> 0");
        assertionMsg(size == 1, size);
        assertion(mesh->vertices().size() == 1);
        ids[0] = 0;
      }
      else {
        preciceDebug("Get read ids");
        utils::DynVector internalPosition(_dimensions);
        utils::DynVector position(_dimensions);
        assertion2(mesh->vertices().size() <= size, mesh->vertices().size(), size);
        for (int i=0; i < size; i++){
          for (int dim=0; dim < _dimensions; dim++){
            position[dim] = positions[i*_dimensions+dim];
          }
          int j=0;
          for (; j < mesh->vertices().size(); j++){
            internalPosition = mesh->vertices()[j].getCoords();
            if (equals(internalPosition, position)){
              ids[i] = j;
              break;
            }
          }
          preciceCheck(j < mesh->vertices().size(), "getReadIDsFromPositions()",
                       "Position " << i << "=" << position << " unknown!");
        }
      }
    }
  }
}

int SolverInterfaceImpl:: getReadNodesSize
(
  int meshID )
{
  preciceTrace1("getReadNodesSize()", meshID);
  if (_clientMode){
    return _requestManager->requestGetReadNodesSize(meshID);
  }
  else {
    MeshContext& context = _accessor->meshContext(meshID);
    if (context.readMappingContext.mapping.get() == NULL){
      return 0;
    }
    else {
      mesh::PtrMesh mesh(context.readMappingContext.localMesh);
      return mesh->vertices().size();
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
    foreach (mesh::Edge& edge, mesh->edges()){
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
    foreach (mesh::Edge& edge, mesh->edges()){
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

void SolverInterfaceImpl:: mapWrittenData
(
  int meshID )
{
  preciceTrace1("mapWrittenData()", meshID);
  if (_clientMode){
    _requestManager->requestMapWrittenData(meshID);
    return;
  }
  impl::MeshContext& context = _accessor->meshContext(meshID);
  impl::MappingContext& mappingContext = context.writeMappingContext;
  if (mappingContext.mapping.get() == NULL){
    preciceWarning("mapWrittenData()", "Mesh \"" << context.mesh->getName()
                   << "\" has no write data to be mapped!");
    return;
  }
  else if (mappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
    preciceWarning("mapWrittenData()", "Mesh \"" << context.mesh->getName()
                   << "\" with incremental write mapping needs not to be mapped!");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    preciceDebug("Compute mapping for mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  foreach (impl::DataContext& context, _accessor->writeDataContexts()){
    if (context.mesh->getID() == meshID){
      int inDataID = context.localData->getID();
      int outDataID = context.data->getID();
      assign(context.data->values()) = 0.0;
      preciceDebug("Map data \"" << context.data->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      mappingContext.mapping->map(inDataID, outDataID);
    }
  }
  mappingContext.hasMappedData = true;
}


void SolverInterfaceImpl:: mapReadData
(
  int meshID )
{
  preciceTrace1 ("mapReadData(int)", meshID);
  if (_clientMode){
    _requestManager->requestMapReadData(meshID);
    return;
  }
  impl::MeshContext& context = _accessor->meshContext(meshID);
  impl::MappingContext& mappingContext = context.readMappingContext;
  if (mappingContext.mapping.get() == NULL){
    preciceWarning("mapReadData()", "Mesh \"" << context.mesh->getName()
                   << "\" has no read data to be mapped!");
    return;
  }
  else if (mappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
    preciceWarning("mapReadData()", "Mesh \"" << context.mesh->getName()
                   << "\" with incremental read mapping needs not to be mapped!");
    return;
  }
  if (not mappingContext.mapping->hasComputedMapping()){
    preciceDebug("Compute read mapping for mesh \"" << context.mesh->getName() << "\"");
    mappingContext.mapping->computeMapping();
  }
  foreach (impl::DataContext& context, _accessor->readDataContexts()){
    if (context.mesh->getID() == meshID){
      int inDataID = context.data->getID();
      int outDataID = context.localData->getID();
      assign(context.localData->values()) = 0.0;
      preciceDebug("Map read data \"" << context.data->getName()
                   << "\" from mesh \"" << context.mesh->getName() << "\"");
      mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.localData->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.localData->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str());
#     endif
    }
  }
  mappingContext.hasMappedData = true;
}

void SolverInterfaceImpl:: writeBlockVectorData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("writeBlockVectorData()", dataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestWriteBlockVectorData(dataID, size, valueIndices, values);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    impl::MappingContext& mapContext = context.mappingContext;
    if (mapContext.mapping.get() == NULL){
      utils::DynVector& valuesInternal = context.data->values();
      for (int i=0; i < size; i++){
        int offsetInternal = valueIndices[i] * _dimensions;
        int offset = i*_dimensions;
        for (int dim=0; dim < _dimensions; dim++){
          assertion2(offset+dim < valuesInternal.size(),
                     offset+dim, valuesInternal.size());
          valuesInternal[offsetInternal + dim] = values[offset + dim];
        }
      }
    }
    else {
      preciceCheck(mapContext.timing != mapping::MappingConfiguration::INCREMENTAL,
                   "writeBlockVectorData()",
                   "Writing block vector data cannot be used with incremental "
                   << "mapping!");
      assertion(context.localData.get() != NULL);
      utils::DynVector& valuesInternal = context.localData->values();
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
}

void SolverInterfaceImpl:: writeVectorData
(
  int           dataID,
  int           valueIndex,
  const double* value )
{
  preciceTrace2 ( "writeVectorData()", dataID, valueIndex );
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
    _requestManager->requestWriteVectorData(dataID, valueIndex, tarch::la::raw(valueCopy));
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    assertion(context.localData.get() != NULL);
    impl::MappingContext& mapContext = context.mappingContext;
    utils::DynVector& values = context.localData->values();
    if ((mapContext.mapping.get() != NULL)
        && (mapContext.timing == mapping::MappingConfiguration::INCREMENTAL)){
      preciceDebug("Map incrementally");
      for (int dim=0; dim < _dimensions; dim++){
        values[dim] = value[dim];
      }
      mapContext.mapping->map(context.localData->getID(), dataID);
    }
    else {
      preciceDebug("Write value directly");
      assertion1(valueIndex >= 0, valueIndex);
      int offset = valueIndex * _dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        values[offset+dim] = value[dim];
      }
    }
  }
}

void SolverInterfaceImpl:: writeBlockScalarData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("writeBlockScalarData()", dataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestWriteBlockScalarData(dataID, size, valueIndices, values);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    impl::MappingContext& mapContext = context.mappingContext;
    if (mapContext.mapping.get() == NULL){
      utils::DynVector& valuesInternal = context.data->values();
      for (int i=0; i < size; i++){
        assertion2(i < valuesInternal.size(), i, valuesInternal.size());
        valuesInternal[valueIndices[i]] = values[i];
      }
    }
    else {
      preciceCheck(mapContext.timing != mapping::MappingConfiguration::INCREMENTAL,
                   "writeBlockScalarData()",
                   "Writing block scalar data cannot be used with incremental "
                   << "mapping!");
      assertion(context.localData.get() != NULL);
      utils::DynVector& valuesInternal = context.localData->values();
      for (int i=0; i < size; i++){
        assertion2(i < valuesInternal.size(), i, valuesInternal.size());
        valuesInternal[valueIndices[i]] = values[i];
      }
    }
  }
}

void SolverInterfaceImpl:: writeScalarData
(
  int    dataID,
  int    valueIndex,
  double value )
{
  preciceTrace3("writeScalarData()", dataID, valueIndex, value );
  preciceCheck(valueIndex >= -1, "writeScalarData()", "Invalid value index ("
               << valueIndex << ") when writing scalar data!");
  if (_clientMode){
    _requestManager->requestWriteScalarData(dataID, valueIndex, value);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    assertion(context.localData.use_count() > 0);
    impl::MappingContext& mapContext = context.mappingContext;
    utils::DynVector& values = context.localData->values();
    bool hasMapping = mapContext.mapping.get() != NULL;
    bool isIncremental = mapContext.timing == mapping::MappingConfiguration::INCREMENTAL;
    if (hasMapping && isIncremental){
      preciceDebug("Map incrementally");
      values[0] = value;
      mapContext.mapping->map(context.localData->getID(), dataID);
    }
    else {
      preciceDebug("Write value directly");
      assertion1(valueIndex >= 0, valueIndex);
      values[valueIndex] = value;
    }
  }
}

void SolverInterfaceImpl:: readBlockVectorData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("readBlockVectorData()", dataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestReadBlockVectorData(dataID, size, valueIndices, values);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    impl::MappingContext& mapContext = context.mappingContext;
    if (mapContext.mapping.get() == NULL){
      utils::DynVector& valuesInternal = context.data->values();
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
    else {
      preciceCheck(mapContext.timing != mapping::MappingConfiguration::INCREMENTAL,
                   "readBlockVectorData()",
                   "Reading block vector data cannot be used with incremental "
                   << "mapping!");
      assertion(context.localData.get() != NULL);
      utils::DynVector& valuesInternal = context.localData->values();
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
}

void SolverInterfaceImpl:: readVectorData
(
  int     dataID,
  int     valueIndex,
  double* value )
{
  preciceTrace2("readVectorData()", dataID, valueIndex);
  preciceCheck(valueIndex >= -1, "readData(vector)", "Invalid value index ( "
               << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadVectorData(dataID, valueIndex, value);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    assertion(context.localData.use_count() > 0);
    utils::DynVector& values = context.localData->values();
    impl::MappingContext& mapContext = context.mappingContext;
    bool hasMapping = mapContext.mapping.get() != NULL;
    bool isIncremental = mapContext.timing == mapping::MappingConfiguration::INCREMENTAL;
    if (hasMapping && isIncremental){
      preciceDebug("Map incrementally");
      assign(context.localData->values()) = 0.0;
      mapContext.mapping->map(dataID, context.localData->getID());
      for (int dim=0; dim < _dimensions; dim++){
        value[dim] = values[dim];
      }
    }
    else {
      preciceDebug("Read (mapped) value directly");
      assertion1 (valueIndex >= 0, valueIndex);
      int offset = valueIndex * _dimensions;
      for (int dim=0; dim < _dimensions; dim++){
        value[dim] = values[offset + dim];
      }
    }
  }
# ifdef Debug
  if (_dimensions == 2) preciceDebug("read value = " << tarch::la::wrap<2>(value));
  if (_dimensions == 3) preciceDebug("read value = " << tarch::la::wrap<3>(value));
# endif
}

void SolverInterfaceImpl:: readBlockScalarData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  preciceTrace2("readBlockScalarData()", dataID, size);
  assertion(valueIndices != NULL);
  assertion(values != NULL);
  if (_clientMode){
    _requestManager->requestReadBlockScalarData(dataID, size, valueIndices, values);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    impl::MappingContext& mapContext = context.mappingContext;
    if (mapContext.mapping.get() == NULL){
      utils::DynVector& valuesInternal = context.data->values();
      for (int i=0; i < size; i++){
        assertion2(valueIndices[i] < valuesInternal.size(),
        		   valueIndices[i], valuesInternal.size());
        values[i] = valuesInternal[valueIndices[i]];
      }
    }
    else {
      preciceCheck(mapContext.timing != mapping::MappingConfiguration::INCREMENTAL,
                   "readBlockScalarData()",
                   "Reading block scalar data cannot be used with incremental "
                   << "mapping!");
      assertion(context.localData.get() != NULL);
      utils::DynVector& valuesInternal = context.localData->values();
      for (int i=0; i < size; i++){
        assertion2(valueIndices[i] < valuesInternal.size(),
        	       valueIndices[i], valuesInternal.size());
        values[i] = valuesInternal[valueIndices[i]];
      }
    }
  }
}

void SolverInterfaceImpl:: readScalarData
(
  int     dataID,
  int     valueIndex,
  double& value )
{
  preciceTrace3("readScalarData()", dataID, valueIndex, value);
  preciceCheck(valueIndex >= -1, "readData(vector)", "Invalid value index ( "
               << valueIndex << " )when reading vector data!");
  if (_clientMode){
    _requestManager->requestReadScalarData(dataID, valueIndex, value);
  }
  else {
    DataContext& context = _accessor->dataContext(dataID);
    assertion(context.localData.use_count() > 0);
    utils::DynVector& values = context.localData->values();
    impl::MappingContext& mapContext = context.mappingContext;
    bool hasMapping = mapContext.mapping.get() != NULL;
    bool isIncremental = mapContext.timing == mapping::MappingConfiguration::INCREMENTAL;
    if (hasMapping && isIncremental){
      preciceDebug("Map incrementally");
      assign(context.localData->values()) = 0.0;
      mapContext.mapping->map(dataID, context.localData->getID());
      value = values[0];
    }
    else {
      preciceDebug("Read (mapped) value directly");
      value = values[valueIndex];
    }
  }
  preciceDebug("Read value = " << value);
}

//void SolverInterfaceImpl:: setExportLocation
//(
//  const std::string& location,
//  int                exportType )
//{
//  preciceTrace2 ( "setExportLocation()", location, exportType );
//  assertion ( not _clientMode ); // TODO implement
//  const ExportContext& context = _accessor->exportContext ();
//  foreach ( const io::PtrExport & exporter, context.exports ) {
//    bool exportAll = exportType == constants::exportAll();
//    bool exportThis = exporter->getType() == exportType;
//    if ( exportAll || exportThis ) {
//      exporter->setLocation ( location );
//    }
//  }
//}

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
  foreach ( const io::ExportContext& context, _accessor->exportContexts() ){
    preciceDebug ( "Export type = " << exportType );
    bool exportAll = exportType == constants::exportAll();
    bool exportThis = context.exporter->getType() == exportType;
    if ( exportAll || exportThis ){
      foreach ( MeshContext& meshContext, _accessor->usedMeshContexts() ){
        std::string name = meshContext.mesh->getName() + "-" + filenameSuffix;
        std::string filename = context.location + name;
        preciceDebug ( "Exporting mesh to file \"" << filename << "\"" );
        context.exporter->doExport ( filename, *meshContext.mesh );
      }
    }
    // Export spacetrees
    if (context.exportSpacetree){
      foreach ( MeshContext& meshContext, _accessor->usedMeshContexts() ){
        std::string name = meshContext.mesh->getName() + "-" + filenameSuffix;
        std::string filename = context.location + name + ".spacetree";
        if ( meshContext.spacetree.get() != NULL ) {
          spacetree::ExportSpacetree exportSpacetree(filename);
          exportSpacetree.doExport ( *meshContext.spacetree );
        }
      }
    }
  }
  // Export neighbors
  //if ( context.plotNeighbors ) {
  //  _exportVTKNeighbors.exportNeighbors ( filenameSuffix + ".neighbors" );
  //}
}

void SolverInterfaceImpl:: integrateData
(
  int     dataID,
  double& integratedValue )
{
  preciceTrace1 ( "integrateData(double)", dataID );
  if (_clientMode){
    _requestManager->requestIntegrateScalarData ( dataID, integratedValue );
  }
  else {
    const DataContext& context = _accessor->dataContext(dataID);
    integratedValue = tarch::la::sum(context.data->values());
  }
}

void SolverInterfaceImpl:: integrateData
(
  int     dataID,
  double* integratedValue )
{
  preciceTrace1("integrateData(Vector)", dataID);
  if (_clientMode){
    _requestManager->requestIntegrateVectorData(dataID, integratedValue);
  }
  else {
    const utils::DynVector& values = _accessor->dataContext(dataID).data->values();
    for (int dim=0; dim < _dimensions; dim++){
      integratedValue[dim] = 0.0;
    }
    for (int i=0; i < values.size(); i += _dimensions){
      for (int dim=0; dim < _dimensions; dim++){
        integratedValue[dim] += values[i+dim];
      }
    }
  }
}

MeshHandle SolverInterfaceImpl:: getMeshHandle
(
  const std::string& meshName )
{
  preciceTrace1("getMeshHandle()", meshName);
  assertion(not _clientMode); // TODO implement
  foreach (MeshContext & context, _accessor->usedMeshContexts()){
    if (context.mesh->getName() == meshName){
      return MeshHandle(context.mesh->content());
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

void SolverInterfaceImpl:: configureCommunications
(
  const com::PtrCommunicationConfiguration& config )
{
  preciceTrace("configureCommunications()");
  typedef com::CommunicationConfiguration::ComTuple ComTuple;
  foreach (ComTuple comTuple, config->communications()){
    std::string comPartner("");
    bool isRequesting = false;
    if (boost::get<1>(comTuple) == _accessorName){
      comPartner = boost::get<2>(comTuple);
      isRequesting = true;
    }
    else if (boost::get<2>(comTuple) == _accessorName){
      comPartner = boost::get<1>(comTuple);
    }
    if (not comPartner.empty()){
      foreach (const impl::PtrParticipant& participant, _participants){
        if (participant->getName() == comPartner){
          if (participant->useServer()){
            comPartner += "Server";
          }
          assertion1(not utils::contained(comPartner, _communications), comPartner);
          assertion(boost::get<0>(comTuple).use_count() > 0);
          Communication com;
          com.communication = boost::get<0>(comTuple);
          com.isRequesting = isRequesting;
          _communications[comPartner] = com;
        }
      }
    }
  }
}

void SolverInterfaceImpl:: configureCommunicatedGeometries
(
  const com::PtrCommunicationConfiguration& comConfig )
{
  preciceTrace ( "configureCommunicatedGeometries()" );
  foreach ( MeshContext& context, _accessor->usedMeshContexts() ) {
    if ( context.provideMesh ) { // Accessor provides geometry
      preciceCheck ( context.receiveMeshFrom.empty(), "configureCommunicatedGeometries()",
                     "Participant \"" << _accessorName << "\" cannot provide "
                     << "and receive mesh " << context.mesh->getName() << "!" );
      utils::DynVector offset ( _dimensions, 0.0 );
      std::string provider ( _accessorName );
      geometry::CommunicatedGeometry* comGeo =
          new geometry::CommunicatedGeometry ( offset, provider, provider );
      bool addedReceiver = false;
      foreach ( PtrParticipant receiver, _participants ){
        foreach ( MeshContext& receiverContext, receiver->usedMeshContexts() ){
          bool doesReceive = receiverContext.receiveMeshFrom == _accessorName;
          doesReceive &= receiverContext.mesh->getName() == context.mesh->getName();
          if ( doesReceive ){
            preciceDebug ( "   ... receiver " << receiver );
            com::PtrCommunication com =
                comConfig->getCommunication ( receiver->getName(), provider );
            comGeo->addReceiver ( receiver->getName(), com );
            addedReceiver = true;
          }
        }
      }
      preciceCheck ( addedReceiver, "configureCommunicatedGeometries()", "No receivers "
                     << " defined for mesh \"" << context.mesh->getName()
                     << "\" provided by participant \"" << _accessorName << "!" );
      preciceCheck ( context.geometry.use_count() == 0, "configureCommunicatedGeometries()",
                     "Participant \"" << _accessorName << "\" cannot provide "
                     << "the geometry of mesh \"" << context.mesh->getName()
                     << " in addition to a defined geometry!" );
      context.geometry = geometry::PtrGeometry ( comGeo );
    }
    else if ( not context.receiveMeshFrom.empty() ) { // Accessor receives geometry
      preciceCheck ( not context.provideMesh, "configureCommunicatedGeometries()",
                     "Participant \"" << _accessorName << "\" cannot provide "
                     << "and receive mesh " << context.mesh->getName() << "!" );
      utils::DynVector offset ( _dimensions, 0.0 );
      std::string receiver ( _accessorName );
      std::string provider ( context.receiveMeshFrom );
      preciceDebug ( "Receiving mesh from " << provider );
      geometry::CommunicatedGeometry * comGeo =
          new geometry::CommunicatedGeometry ( offset, receiver, provider );
      com::PtrCommunication com = comConfig->getCommunication ( receiver, provider );
      comGeo->addReceiver ( receiver, com );
      preciceCheck ( context.geometry.use_count() == 0, "configureCommunicatedGeometries()",
                     "Participant \"" << _accessorName << "\" cannot receive "
                     << "the geometry of mesh \"" << context.mesh->getName()
                     << " in addition to a defined geometry!" );
      context.geometry = geometry::PtrGeometry ( comGeo );
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
        geometry::ImportGeometry::VRML_1_FILE, true);
    geometry = geometry::PtrGeometry ( importGeo );
    bool importLocal = meshContext.writeMappingContext.localMesh.get() != NULL;
    importLocal &= meshContext.writeMappingContext.timing != mapping::MappingConfiguration::INCREMENTAL;
    if (importLocal){
      fileName = "precice_checkpoint_" + _accessorName + "_" + meshName + "_localwrite";
      geometry::ImportGeometry importGeo(utils::DynVector(_dimensions, 0.0),
          fileName, geometry::ImportGeometry::VRML_1_FILE, true);
      importGeo.create(*meshContext.writeMappingContext.localMesh);
    }
    importLocal = meshContext.readMappingContext.localMesh.get() != NULL;
    importLocal &= meshContext.readMappingContext.timing != mapping::MappingConfiguration::INCREMENTAL;
    if (importLocal){
      fileName = "precice_checkpoint_" + _accessorName + "_" + meshName + "_localread";
      geometry::ImportGeometry importGeo(utils::DynVector(_dimensions, 0.0),
          fileName, geometry::ImportGeometry::VRML_1_FILE, true);
      importGeo.create(*meshContext.readMappingContext.localMesh);
    }
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

  // Create default vertex for temporary participant meshes
  if (meshContext.writeMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
    mesh::PtrMesh& localMesh = meshContext.writeMappingContext.localMesh;
    assertion(localMesh != meshContext.mesh);
    assertion(localMesh->vertices().size() == 0);
    localMesh->createVertex(utils::DynVector(_dimensions,0.0));
    localMesh->allocateDataValues();
  }
  if (meshContext.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL){
    mesh::PtrMesh & localMesh = meshContext.readMappingContext.localMesh;
    assertion(localMesh != meshContext.mesh);
    assertion(localMesh->vertices().size() == 0);
    localMesh->createVertex(utils::DynVector(_dimensions,0.0));
    localMesh->allocateDataValues();
  }
}

void SolverInterfaceImpl:: mapWrittenData()
{
  preciceTrace("mapWrittenData()");
  using namespace mapping;
  MappingConfiguration::Timing timing;
  // Compute mappings
  foreach (impl::MeshContext& context, _accessor->usedMeshContexts()){
    bool hasMapping = context.writeMappingContext.mapping.get() != NULL;
    timing = context.writeMappingContext.timing;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    if (hasMapping && rightTime){
      bool hasComputed = context.writeMappingContext.mapping->hasComputedMapping();
      if (not hasComputed){
        preciceDebug("Compute write mapping for mesh \"" << context.mesh->getName() << "\"");
        context.writeMappingContext.mapping->computeMapping();
      }
    }
  }
  // Map data
  foreach (impl::DataContext& context, _accessor->writeDataContexts()){
    timing = context.mappingContext.timing;
    bool hasMapping = context.mappingContext.mapping.get() != NULL;
    bool rightTime = timing == MappingConfiguration::ON_ADVANCE;
    rightTime |= timing == MappingConfiguration::INITIAL;
    bool hasMapped = context.mappingContext.hasMappedData;
    if (hasMapping && rightTime && (not hasMapped)){
      int inDataID = context.localData->getID();
      int outDataID = context.data->getID();
      preciceDebug("Map data \"" << context.data->getName()
                   << "\" to mesh \"" << context.mesh->getName() << "\"");
      assign(context.data->values()) = 0.0;
      context.mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.data->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.data->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str() );
#     endif
    }
  }

  // Clear non-stationary, non-incremental mappings
  foreach (impl::MeshContext& context, _accessor->usedMeshContexts()){
    bool hasMapping = context.writeMappingContext.mapping.use_count() > 0;
    bool isIncremental = context.writeMappingContext.timing
                         == MappingConfiguration::INCREMENTAL;
    bool isStationary = context.writeMappingContext.timing
                        == MappingConfiguration::INITIAL;
    if (hasMapping && (not isIncremental) && (not isStationary)){
        context.writeMappingContext.mapping->clear();
    }
    context.writeMappingContext.hasMappedData = false;
  }
}

void SolverInterfaceImpl:: mapReadData()
{
  preciceTrace("mapReadData()");
  mapping::MappingConfiguration::Timing timing;
  // Compute mappings
  foreach (impl::MeshContext& context, _accessor->usedMeshContexts()){
    timing = context.readMappingContext.timing;
    bool mapNow = timing == mapping::MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == mapping::MappingConfiguration::INITIAL;
    bool hasMapping = context.readMappingContext.mapping.use_count() > 0;
    bool isIncremental = context.readMappingContext.timing
                         == mapping::MappingConfiguration::INCREMENTAL;
    if (mapNow && hasMapping && (not isIncremental)){
      // Compute mapping only when the mapping has not been computed. For timing
      // initial, this is the case only once.
      if (not context.readMappingContext.mapping->hasComputedMapping()){
        preciceDebug("Compute read mapping for mesh \"" << context.mesh->getName() << "\"");
        context.readMappingContext.mapping->computeMapping();
      }
    }
  }
  // Map data
  foreach (impl::DataContext& context, _accessor->readDataContexts()){
    timing = context.mappingContext.timing;
    bool mapNow = timing == mapping::MappingConfiguration::ON_ADVANCE;
    mapNow |= timing == mapping::MappingConfiguration::INITIAL;
    bool hasMapping = context.mappingContext.mapping.get() != NULL;
    bool isIncremental = context.mappingContext.timing == mapping::MappingConfiguration::INCREMENTAL;
    bool hasMapped = context.mappingContext.hasMappedData;
    if (mapNow && hasMapping && (not isIncremental) && (not hasMapped)){
      int inDataID = context.data->getID();
      int outDataID = context.localData->getID();
      assign(context.localData->values()) = 0.0;
      preciceDebug("Map read data \"" << context.data->getName()
                   << "\" from mesh \"" << context.mesh->getName() << "\"");
      context.mappingContext.mapping->map(inDataID, outDataID);
#     ifdef Debug
      int max = context.localData->values().size();
      std::ostringstream stream;
      for (int i=0; (i < max) && (i < 10); i++){
        stream << context.localData->values()[i] << " ";
      }
      preciceDebug("First mapped values = " << stream.str());
#     endif
    }
  }

  // Clear non-initial, non-incremental mappings
  foreach (impl::MeshContext& context, _accessor->usedMeshContexts()){
    bool hasMapping = context.readMappingContext.mapping.use_count() > 0;
    if (hasMapping){
      bool isIncremental = context.readMappingContext.timing == mapping::MappingConfiguration::INCREMENTAL;
      bool isStationary = context.readMappingContext.timing
                          == mapping::MappingConfiguration::INITIAL;
      if ((not isIncremental) && (not isStationary)){
        context.readMappingContext.mapping->clear();
      }
      context.readMappingContext.hasMappedData = false;
    }
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
  foreach (action::PtrAction& action, _accessor->actions()){
    if (timings.find(action->getTiming()) != timings.end()){
      action->performAction(time, dt, partFullDt, fullDt);
    }
  }
}

void SolverInterfaceImpl:: handleExports()
{
  preciceTrace("handleExports()");
  assertion(not _clientMode);
  foreach (const io::ExportContext& context, _accessor->exportContexts()){
    if (_couplingScheme->isCouplingTimestepComplete() || context.everyIteration){
      if (context.timestepInterval != -1){
        if (_couplingScheme->getTimesteps() % context.timestepInterval == 0){
          if (context.everyIteration){
            std::ostringstream everySuffix;
            everySuffix << _accessorName << ".it" << _numberAdvanceCalls;
            exportMesh(everySuffix.str());
          }
          std::ostringstream suffix;
          suffix << _accessorName << ".dt" << _couplingScheme->getTimesteps();
          exportMesh(suffix.str());
          if (context.triggerSolverPlot){
            _couplingScheme->requireAction(constants::actionPlotOutput());
          }
        }
      }
    }
  }

  if (_couplingScheme->isCouplingTimestepComplete()){
    // Export watch point data
    foreach (PtrWatchPoint watchPoint, _accessor->watchPoints()){
      watchPoint->exportPointData(_couplingScheme->getTime());
    }

    // Checkpointing
    int checkpointingInterval = _couplingScheme->getCheckpointTimestepInterval();
    int timestep = _couplingScheme->getTimesteps();
    if ((checkpointingInterval != -1) && (timestep % checkpointingInterval == 0)){
      preciceDebug("Set require checkpoint");
      _couplingScheme->requireAction(constants::actionWriteSimulationCheckpoint());
      foreach (const MeshContext& meshContext, _accessor->usedMeshContexts()){
        io::ExportVRML exportVRML(false);
        std::string filename("precice_checkpoint_" + _accessorName
                             + "_" + meshContext.mesh->getName());
        exportVRML.doExportCheckpoint(filename, *meshContext.mesh);
        if (meshContext.writeMappingContext.localMesh.get() != NULL){
          filename = "precice_checkpoint_" + _accessorName + "_" +
                     meshContext.mesh->getName() + "_" + "localwrite";
          exportVRML.doExportCheckpoint(filename, *meshContext.writeMappingContext.localMesh);
        }
        if (meshContext.readMappingContext.localMesh.get() != NULL){
          filename = "precice_checkpoint_" + _accessorName + "_" +
                     meshContext.mesh->getName() + "_" + "localread";
          exportVRML.doExportCheckpoint(filename, *meshContext.readMappingContext.localMesh);
        }
      }
      io::SimulationStateIO exportState(_checkpointFileName + "_simstate.txt");
      exportState.writeState(_couplingScheme->getTime(), timestep, _numberAdvanceCalls);
      //io::TXTWriter exportCouplingSchemeState(_checkpointFileName + "_cplscheme.txt");
      _couplingScheme->exportState(_checkpointFileName);
    }
  }
}

void SolverInterfaceImpl:: resetWrittenData()
{
  preciceTrace("resetWrittenData()");
  foreach (DataContext& context, _accessor->writeDataContexts()){
    assign(context.data->values()) = 0.0;
    if (context.localData != context.data){
      assign(context.localData->values()) = 0.0;
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
  foreach ( const PtrParticipant& participant, partConfig->getParticipants() ) {
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
      const MeshContext& context = _accessor->usedMeshContexts()[i];
      if (context.spacetree.get() == NULL){
        markedMeshContexts[i] = markedQueryDirectly();
      }
      else if (context.mesh->getID() == context.spacetree->meshes().front()->getID()){
        markedMeshContexts[i] = markedQuerySpacetree();
      }
      else {
        markedMeshContexts[i] = markedSkip();
      }
    }
  }
  else {
    for (int i=0; i < (int)markedMeshContexts.size(); i++){
      const MeshContext& context = _accessor->usedMeshContexts()[i];
      if (utils::contained(context.mesh->getID(), meshIDs)){
        if (context.spacetree.get() == NULL){
          markedMeshContexts[i] = markedQueryDirectly();
        }
        else {
          bool allSpacetreeMeshesAreInquired = true;
          foreach (const mesh::PtrMesh& mesh, context.spacetree->meshes()){
            if (not utils::contained(mesh->getID(), meshIDs)){
              allSpacetreeMeshesAreInquired = false;
              break;
            }
          }
          if (allSpacetreeMeshesAreInquired){
            bool isFirst = context.mesh->getID()
                           == context.spacetree->meshes().front()->getID();
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
  com::PtrCommunication com = _accessor->getClientServerCommunication();
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

}} // namespace precice, impl

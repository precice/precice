// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Participant.hpp"
#include "DataContext.hpp"
#include "MeshContext.hpp"
#include "MappingContext.hpp"
#include "action/Action.hpp"
#include "WatchPoint.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "spacetree/Spacetree.hpp"

namespace precice {
namespace impl {

tarch::logging::Log Participant:: _log ( "precice::impl::Participant" );

int Participant:: _participantsSize = 0;

void Participant:: resetParticipantCount()
{
  _participantsSize = 0;
}

Participant:: Participant
(
  const std::string&          name,
  mesh::PtrMeshConfiguration& meshConfig )
:
  _name ( name ),
  _id ( _participantsSize ),
  _watchPoints (),
  _exportContexts(),
  _actions (),
  _meshContexts ( meshConfig->meshes().size(), NULL ),
  _readMappingContexts(),
  _writeMappingContexts(),
  _usedMeshContexts (),
  _dataContexts ( meshConfig->getDataConfiguration()->data().size()*meshConfig->meshes().size(), NULL ),
  _writeDataContexts (),
  _readDataContexts (),
  _clientServerCommunication (),
  _masterSlaveCommunication ()
{
  _participantsSize ++;
}

Participant:: ~Participant()
{
  for(MeshContext* context : _usedMeshContexts){
    delete context;
  }
  _usedMeshContexts.clear();
  _readDataContexts.deleteElements();
  _writeDataContexts.deleteElements();
  _readMappingContexts.deleteElements();
  _writeMappingContexts.deleteElements();
  _participantsSize--;
}

const std::string& Participant:: getName() const
{
  return _name;
}

int Participant:: getID() const
{
  return _id;
}

void Participant:: addWatchPoint
(
  const PtrWatchPoint& watchPoint )
{
  _watchPoints.push_back ( watchPoint );
}

std::vector<PtrWatchPoint>& Participant:: watchPoints()
{
  return _watchPoints;
}

void Participant:: useMesh
(
  const mesh::PtrMesh&                   mesh,
  const geometry::PtrGeometry&   geometry,
  const spacetree::PtrSpacetree& spacetree,
  const utils::DynVector&                localOffset,
  bool                                   remote,
  const std::string&                     fromParticipant,
  double                                 safetyFactor,
  bool                                   provideMesh )
{
  preciceTrace3 ( "useMesh()", _name,  mesh->getName(), mesh->getID() );
  checkDuplicatedUse(mesh);
  assertion ( mesh->getID() < (int)_meshContexts.size() );
  MeshContext* context = new MeshContext(mesh->getDimensions());
  context->mesh = mesh;
  context->geometry = geometry;
  context->spacetree = spacetree;
  context->localOffset = localOffset;
  assertion2 ( mesh->getDimensions() == context->localOffset.size(),
               mesh->getDimensions(), context->localOffset.size() );
  context->receiveMeshFrom = fromParticipant;
  context->safetyFactor = safetyFactor;
  context->provideMesh = provideMesh;

//  if ( spacetree.use_count() > 0 ) {
//    spacetree->setCenter ( spacetree->getCenter() + localOffset );
//  }
//  if ( provideMesh ) {
//    context->meshRequirement = mapping::Mapping::FULL;
//  }

  _meshContexts[mesh->getID()] = context;

  _usedMeshContexts.push_back ( context );

  preciceCheck ( fromParticipant.empty() || (! provideMesh), "useMesh()",
                 "Participant " << _name << " cannot receive and provide mesh "
                 << mesh->getName() << " at the same time!" );
}

void Participant:: addWriteData
(
  const mesh::PtrData& data,
  const mesh::PtrMesh& mesh )
{
  checkDuplicatedData ( data );
  assertion ( data->getID() < (int)_dataContexts.size() );
  DataContext* context = new DataContext ();
  context->fromData = data;
  context->mesh = mesh;
  // will be overwritten later if a mapping exists
  context->toData = context->fromData;
  _dataContexts[data->getID()] = context;
  _writeDataContexts.push_back ( context );
}

void Participant:: addReadData
(
  const mesh::PtrData& data,
  const mesh::PtrMesh& mesh )
{
  checkDuplicatedData ( data );
  assertion ( data->getID() < (int)_dataContexts.size() );
  DataContext* context = new DataContext ();
  context->toData = data;
  context->mesh = mesh;
  // will be overwritten later if a mapping exists
  context->fromData = context->toData;
  _dataContexts[data->getID()] = context;
  _readDataContexts.push_back ( context );
}

void Participant::addReadMappingContext
(
  MappingContext* mappingContext)
{
  _readMappingContexts.push_back(mappingContext);
}

void Participant::addWriteMappingContext
(
  MappingContext* mappingContext)
{
  _writeMappingContexts.push_back(mappingContext);
}

const utils::ptr_vector<MappingContext>& Participant::readMappingContexts() const
{
  return _readMappingContexts;
}

const utils::ptr_vector<MappingContext>& Participant::writeMappingContexts() const
{
  return _writeMappingContexts;
}

const DataContext& Participant:: dataContext
(
  int dataID ) const
{
  assertion ( (dataID >= 0) && (dataID < (int)_dataContexts.size()) );
  assertion ( _dataContexts[dataID] != NULL );
  return *_dataContexts[dataID];
}

DataContext& Participant:: dataContext
(
  int dataID )
{
  preciceTrace2 ( "dataContext(id)", dataID, _dataContexts.size() );
  assertion ( (dataID >= 0) && (dataID < (int)_dataContexts.size()) );
  assertion ( _dataContexts[dataID] != NULL );
  return *_dataContexts[dataID];
}

const utils::ptr_vector<DataContext>& Participant:: writeDataContexts() const
{
  return _writeDataContexts;
}

utils::ptr_vector<DataContext>& Participant:: writeDataContexts()
{
  return _writeDataContexts;
}

const utils::ptr_vector<DataContext>& Participant:: readDataContexts() const
{
  return _readDataContexts;
}

utils::ptr_vector<DataContext>& Participant:: readDataContexts()
{
  return _readDataContexts;
}

bool Participant:: isMeshUsed
(
  int meshID ) const
{
  assertion ( (meshID >= 0) && (meshID < (int)_meshContexts.size()) );
  return _meshContexts[meshID] != NULL;
}

bool Participant:: isDataUsed
(
  int dataID ) const
{
  assertion ( (dataID >= 0) && (dataID < (int)_dataContexts.size()) );
  return _dataContexts[dataID] != NULL;
}

const MeshContext& Participant:: meshContext
(
  int meshID ) const
{
  assertion((meshID >= 0) && (meshID < (int)_meshContexts.size()));
  assertion(_meshContexts[meshID] != NULL);
  return *_meshContexts[meshID];
}

MeshContext& Participant:: meshContext
(
  int meshID )
{
  preciceTrace2("meshContext()", meshID, _meshContexts.size());
  assertion2((meshID >= 0) && (meshID < (int)_meshContexts.size()),
             meshID, _meshContexts.size());
  assertion(_meshContexts[meshID] != NULL);
  return *_meshContexts[meshID];
}

const std::vector<MeshContext*>& Participant:: usedMeshContexts() const
{
  return _usedMeshContexts;
}

std::vector<MeshContext*>& Participant:: usedMeshContexts()
{
  return _usedMeshContexts;
}

void Participant:: addAction
(
  const action::PtrAction& action )
{
  _actions.push_back ( action );
}

std::vector<action::PtrAction>& Participant:: actions()
{
  return _actions;
}

const std::vector<action::PtrAction>& Participant:: actions() const
{
  return _actions;
}

void Participant:: addExportContext
(
  const io::ExportContext& exportContext )
{
  _exportContexts.push_back(exportContext);
}

const std::vector<io::ExportContext>& Participant:: exportContexts() const
{
  return _exportContexts;
}

void Participant:: checkDuplicatedUse
(
  const mesh::PtrMesh& mesh )
{
  assertion ( (int)_meshContexts.size() > mesh->getID() );
  preciceCheck ( _meshContexts[mesh->getID()] == NULL, "checkDuplicateUse()",
                 "Mesh \"" << mesh->getName() << " cannot be used twice by "
                 << "participant " << _name << "!" );
}

void Participant:: checkDuplicatedData
(
  const mesh::PtrData& data )
{
  preciceTrace2 ( "checkDuplicatedData()", data->getID(), _dataContexts.size() );
  assertion2 ( data->getID() < (int)_dataContexts.size(), data->getID(), _dataContexts.size() );
  preciceCheck ( _dataContexts[data->getID()] == NULL, "checkDuplicatedData()",
                 "Participant \"" << _name << "\" can read/write data \""
                 << data->getName() << "\" only once!" );
}

bool Participant:: useServer()
{
  return _clientServerCommunication.use_count() > 0;
}

void Participant:: setClientServerCommunication
(
  com::Communication::SharedPointer communication )
{
  assertion ( communication.use_count() > 0 );
  _clientServerCommunication = communication;
}

com::Communication::SharedPointer Participant:: getClientServerCommunication() const
{
  return _clientServerCommunication;
}

bool Participant:: useMaster()
{
  return _masterSlaveCommunication.use_count() > 0;
}

void Participant:: setMasterSlaveCommunication
(
  com::Communication::SharedPointer communication )
{
  assertion ( communication.use_count() > 0 );
  _masterSlaveCommunication = communication;
}

com::Communication::SharedPointer Participant:: getMasterSlaveCommunication() const
{
  return _masterSlaveCommunication;
}


}} // namespace precice, impl

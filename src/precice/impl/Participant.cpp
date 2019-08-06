#include "Participant.hpp"
#include "DataContext.hpp"
#include "MeshContext.hpp"
#include "MappingContext.hpp"
#include "action/Action.hpp"
#include "WatchPoint.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include <utility>
#include <algorithm>

namespace precice {
namespace impl {

int Participant:: _participantsSize = 0;

void Participant:: resetParticipantCount()
{
  _participantsSize = 0;
}

Participant:: Participant
(
  std::string                 name,
  mesh::PtrMeshConfiguration& meshConfig )
:
  _name (std::move(name)),
  _id ( _participantsSize ),
  _meshContexts ( meshConfig->meshes().size(), nullptr ),
  _dataContexts ( meshConfig->getDataConfiguration()->data().size()*meshConfig->meshes().size(), nullptr )
{
  _participantsSize++;
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
  const mesh::PtrMesh&                          mesh,
  const Eigen::VectorXd&                        localOffset,
  bool                                          remote,
  const std::string&                            fromParticipant,
  double                                        safetyFactor,
  bool                                          provideMesh,
  partition::ReceivedPartition::GeometricFilter geoFilter)
{
  P_TRACE(_name,  mesh->getName(), mesh->getID() );
  checkDuplicatedUse(mesh);
  P_ASSERT( mesh->getID() < (int)_meshContexts.size() );
  auto context = new MeshContext(mesh->getDimensions());
  context->mesh = mesh;
  context->localOffset = localOffset;
  P_ASSERT( mesh->getDimensions() == context->localOffset.size(),
               mesh->getDimensions(), context->localOffset.size() );
  context->receiveMeshFrom = fromParticipant;
  context->safetyFactor = safetyFactor;
  context->provideMesh = provideMesh;
  context->geoFilter = geoFilter;

  _meshContexts[mesh->getID()] = context;

  _usedMeshContexts.push_back ( context );

  P_CHECK( fromParticipant.empty() || (! provideMesh),
         "Participant " << _name << " cannot receive and provide mesh "
         << mesh->getName() << " at the same time!" );
}

void Participant:: addWriteData
(
  const mesh::PtrData& data,
  const mesh::PtrMesh& mesh )
{
  checkDuplicatedData ( data );
  P_ASSERT( data->getID() < (int)_dataContexts.size() );
  auto context = new DataContext ();
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
  P_ASSERT( data->getID() < (int)_dataContexts.size() );
  auto context = new DataContext ();
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
  P_ASSERT( (dataID >= 0) && (dataID < (int)_dataContexts.size()) );
  P_ASSERT( _dataContexts[dataID] != nullptr );
  return *_dataContexts[dataID];
}

DataContext& Participant:: dataContext
(
  int dataID )
{
  P_TRACE(dataID, _dataContexts.size() );
  P_ASSERT( (dataID >= 0) && (dataID < (int)_dataContexts.size()) );
  P_ASSERT( _dataContexts[dataID] != nullptr );
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
  P_ASSERT( (meshID >= 0) && (meshID < (int)_meshContexts.size()) );
  return _meshContexts[meshID] != nullptr;
}

bool Participant:: isDataUsed
(
  int dataID ) const
{
  P_ASSERT( (dataID >= 0) && (dataID < (int)_dataContexts.size()), dataID, (int)_dataContexts.size() );
  return _dataContexts[dataID] != nullptr;
}

bool Participant:: isDataRead
(
  int dataID ) const
{
    return std::any_of(_readDataContexts.begin(), _readDataContexts.end(), [dataID](const DataContext& context) {
            return context.toData->getID() == dataID;
            });
}

bool Participant:: isDataWrite
(
  int dataID ) const
{
    return std::any_of(_writeDataContexts.begin(), _writeDataContexts.end(), [dataID](const DataContext& context) {
            return context.fromData->getID() == dataID;
            });
}

const MeshContext& Participant:: meshContext
(
  int meshID ) const
{
  P_ASSERT((meshID >= 0) && (meshID < (int)_meshContexts.size()));
  P_ASSERT(_meshContexts[meshID] != nullptr);
  return *_meshContexts[meshID];
}

MeshContext& Participant:: meshContext
(
  int meshID )
{
  P_TRACE(meshID, _meshContexts.size());
  P_ASSERT((meshID >= 0) && (meshID < (int)_meshContexts.size()),
             meshID, _meshContexts.size());
  P_ASSERT(_meshContexts[meshID] != nullptr);
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
  auto& context = meshContext(action->getMesh()->getID());
  context.require(action->getMeshRequirement());
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
  P_ASSERT( (int)_meshContexts.size() > mesh->getID() );
  P_CHECK( _meshContexts[mesh->getID()] == nullptr,
         "Mesh \"" << mesh->getName() << " cannot be used twice by "
         << "participant " << _name << "!" );
}

void Participant:: checkDuplicatedData
(
  const mesh::PtrData& data )
{
  P_TRACE(data->getID(), _dataContexts.size() );
  P_ASSERT( data->getID() < (int)_dataContexts.size(), data->getID(), _dataContexts.size() );
  P_CHECK( _dataContexts[data->getID()] == nullptr,
          "Participant \"" << _name << "\" can read/write data \""
          << data->getName() << "\" only once!" );
}

bool Participant:: useServer()
{
  return static_cast<bool>(_clientServerCommunication);
}

void Participant:: setClientServerCommunication
(
  com::PtrCommunication communication )
{
  P_ASSERT( communication );
  _clientServerCommunication = std::move(communication);
}

com::PtrCommunication Participant:: getClientServerCommunication() const
{
  return _clientServerCommunication;
}

bool Participant:: useMaster()
{
  return _useMaster;
}

void Participant:: setUseMaster(bool useMaster)
{
  _useMaster = useMaster;
}


}} // namespace precice, impl

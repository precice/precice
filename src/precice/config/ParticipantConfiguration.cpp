// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParticipantConfiguration.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/MappingContext.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "spacetree/config/SpacetreeConfiguration.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Dimensions.hpp"

namespace precice {
namespace config {

tarch::logging::Log ParticipantConfiguration::
    _log("precice::config::ParticipantConfiguration");

//const std::string& ParticipantConfiguration:: getTag()
//{
//  static std::string tag("participant");
//  return tag;
//}

ParticipantConfiguration:: ParticipantConfiguration
(
  utils::XMLTag&                              parent,
  const mesh::PtrMeshConfiguration&           meshConfiguration,
  const geometry::PtrGeometryConfiguration&   geometryConfiguration,
  const spacetree::PtrSpacetreeConfiguration& spacetreeConfiguration )
:
  TAG("participant"),
  TAG_WRITE("write-data"),
  TAG_READ("read-data"),
  TAG_DATA_ACTION("data-action"),
  TAG_USE_MESH("use-mesh"),
  TAG_WATCH_POINT("watch-point"),
  TAG_SERVER("server"),
  TAG_MASTER("master"),
  ATTR_NAME("name"),
  ATTR_SOURCE_DATA("source-data"),
  ATTR_TARGET_DATA("target-data"),
  ATTR_TIMING("timing"),
  ATTR_LOCAL_OFFSET("offset"),
  ATTR_ACTION_TYPE("type"),
  ATTR_FROM("from"),
  ATTR_SAFETY_FACTOR("safety-factor"),
  ATTR_PROVIDE("provide"),
  ATTR_MESH("mesh"),
  ATTR_COORDINATE("coordinate"),
  ATTR_COMMUNICATION("communication"),
  ATTR_CONTEXT("context"),
  ATTR_NETWORK("network"),
  ATTR_EXCHANGE_DIRECTORY("exchange-directory"),
  _dimensions(0),
  _meshConfig(meshConfiguration),
  _geometryConfig(geometryConfiguration),
  _spacetreeConfig(spacetreeConfiguration),
  _mappingConfig(),
  _actionConfig(),
  //_isValid(false),
  _participants(),
  _watchPointConfigs()
{
  //assertion1 ( (_dimensions == 2) || (_dimensions == 3), _dimensions );
  assertion(_meshConfig.use_count() > 0);
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Represents one solver using preCICE. In a coupled simulation, two ";
  doc += "participants have to be defined, while in geometry mode (see tag ";
  doc += "<solver-interface>) only one participant is necessary.";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrName(ATTR_NAME);
  doc = "Name of the participant. Has to match the name given on construction ";
  doc += "of the precice::SolverInterface object used by the participant.";
  attrName.setDocumentation(doc);
  tag.addAttribute(attrName);

  XMLTag tagWriteData(*this, TAG_WRITE, XMLTag::OCCUR_ARBITRARY);
  doc = "Sets data to be written by the participant to preCICE. ";
  doc += "Data is defined by using the <data> tag.";
  tagWriteData.setDocumentation(doc);
  XMLTag tagReadData(*this, TAG_READ, XMLTag::OCCUR_ARBITRARY);
  doc = "Sets data to be read by the participant from preCICE. ";
  doc += "Data is defined by using the <data> tag.";
  tagReadData.setDocumentation(doc);
  XMLAttribute<std::string> attrDataName(ATTR_NAME);
  doc = "Name of the data.";
  attrDataName.setDocumentation(doc);
  tagWriteData.addAttribute(attrDataName);
  tagReadData.addAttribute(attrDataName);
  XMLAttribute<std::string> attrMesh(ATTR_MESH);
  doc = "Mesh the data belongs to. If data should be read/written to several ";
  doc += "meshes, this has to be specified separately for each mesh.";
  attrMesh.setDocumentation(doc);
  tagWriteData.addAttribute(attrMesh);
  tagReadData.addAttribute(attrMesh);
  tag.addSubtag(tagWriteData);
  tag.addSubtag(tagReadData);

  _mappingConfig = mapping::PtrMappingConfiguration(
                   new mapping::MappingConfiguration(tag, _meshConfig));

  _actionConfig = action::PtrActionConfiguration(
                  new action::ActionConfiguration(tag, _meshConfig));

  _exportConfig = io::PtrExportConfiguration(new io::ExportConfiguration(tag));

  XMLTag tagWatchPoint(*this, TAG_WATCH_POINT, XMLTag::OCCUR_ARBITRARY);
  doc = "A watch point can be used to follow the transient changes of data ";
  doc += "and mesh vertex coordinates at a given point";
  tagWatchPoint.setDocumentation(doc);
  doc = "Name of the watch point. Is taken in combination with the participant ";
  doc += "name to construct the filename the watch point data is written to.";
  attrName.setDocumentation(doc);
  tagWatchPoint.addAttribute(attrName);
  doc = "Mesh to be watched.";
  attrMesh.setDocumentation(doc);
  tagWatchPoint.addAttribute(attrMesh);
  XMLAttribute<utils::DynVector> attrCoordinate(ATTR_COORDINATE);
  doc = "The coordinates of the watch point. I the watch point is not put exactly ";
  doc += "on the mesh to observe, the closest projection of the point onto the ";
  doc += "mesh is considered instead, and values/coordinates are interpolated ";
  doc += "linearly to that point.";
  attrCoordinate.setDocumentation(doc);
  tagWatchPoint.addAttribute(attrCoordinate);
  tag.addSubtag(tagWatchPoint);

  XMLTag tagUseMesh(*this, TAG_USE_MESH, XMLTag::OCCUR_ARBITRARY);
  doc = "Makes a mesh (see tag <mesh> available to a participant.";
  tagUseMesh.setDocumentation(doc);
  attrName.setDocumentation("Name of the mesh.");
  tagUseMesh.addAttribute(attrName);
  XMLAttribute<utils::DynVector> attrLocalOffset(ATTR_LOCAL_OFFSET);
  doc = "The mesh can have an offset only applied for the local participant. ";
  doc += "Vector-valued example: '1.0; 0.0; 0.0'";
  attrLocalOffset.setDocumentation(doc);
  attrLocalOffset.setDefaultValue(utils::DynVector(3, 0.0));
  tagUseMesh.addAttribute(attrLocalOffset);

  XMLAttribute<std::string> attrFrom(ATTR_FROM);
  doc = "A mesh might not be constructed by a geometry (see tags <geometry:...>), ";
  doc += "but by a solver directly. If a solver created mesh should be used by ";
  doc += "another solver, this attribute has to specify the creating participant's";
  doc += " name. The creator has to use the attribute \"provide\" to signal he is ";
  doc += "providing the mesh geometry.";
  attrFrom.setDocumentation(doc);
  attrFrom.setDefaultValue("");
  tagUseMesh.addAttribute(attrFrom);
  XMLAttribute<double> attrSafetyFactor(ATTR_SAFETY_FACTOR);
  doc = "If a mesh is receiced from another partipant (see tag <from>), it needs to ";
  doc += "decomposed on this participant (in master mode only). To speed up this process, ";
  doc += "the receiving master pre-filters every mesh by bounding box information from the ";
  doc += "local mesh. This safety factor defines how by which factor this local information ";
  doc += "increased. An example: 0.1 means that the bounding box is 110% of its original size.";
  attrSafetyFactor.setDocumentation(doc);
  attrSafetyFactor.setDefaultValue(0.1);
  tagUseMesh.addAttribute(attrSafetyFactor);
  XMLAttribute<bool> attrProvide(ATTR_PROVIDE);
  doc = "A mesh might not be constructed by a geometry (see tags<geometry:...>), ";
  doc += "but by a solver directly. If this attribute is set to \"on\", the ";
  doc += "participant has to create the mesh geometry before initializing preCICE.";
  attrProvide.setDocumentation(doc);
  attrProvide.setDefaultValue(false);
  tagUseMesh.addAttribute(attrProvide);
  tag.addSubtag(tagUseMesh);

  std::list<XMLTag> serverTags;
  XMLTag::Occurrence serverOcc = XMLTag::OCCUR_NOT_OR_ONCE;
  {
    XMLTag tagServer(*this, "sockets", serverOcc, TAG_SERVER);
    doc = "When a solver runs in parallel, it has to use preCICE in form of a ";
    doc += "separatly running server. This is enabled by this tag. ";
    doc += "The communication between participant and server is done by sockets.";
    tagServer.setDocumentation(doc);

    XMLAttribute<int> attrPort("port");
    doc = "Port number (16-bit unsigned integer) to be used for socket ";
    doc += "communiation. The default is \"0\", what means that OS will ";
    doc += "dynamically search for a free port (if at least one exists) and ";
    doc += "bind it automatically.";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(0);
    tagServer.addAttribute(attrPort);

    XMLAttribute<std::string> attrNetwork(ATTR_NETWORK);
    doc = "Network name to be used for socket communiation. ";
    doc += "Default is \"lo\", i.e., the local host loopback.";
    attrNetwork.setDocumentation(doc);
    attrNetwork.setDefaultValue("lo");
    tagServer.addAttribute(attrNetwork);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where connection information is exchanged. By default, the ";
    doc += "directory of startup is chosen, and both solvers have to be started ";
    doc += "in the same directory.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tagServer.addAttribute(attrExchangeDirectory);

    serverTags.push_back(tagServer);
  }
  {
    XMLTag tagServer(*this, "mpi", serverOcc, TAG_SERVER);
    doc = "When a solver runs in parallel, it has to use preCICE in form of a ";
    doc += "separatly running server. This is enabled by this tag. ";
    doc += "The communication between participant and server is done by mpi ";
    doc += "with startup in separated communication spaces.";
    tagServer.setDocumentation(doc);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where connection information is exchanged. By default, the ";
    doc += "directory of startup is chosen, and both solvers have to be started ";
    doc += "in the same directory.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tagServer.addAttribute(attrExchangeDirectory);

    serverTags.push_back(tagServer);
  }
  {
    XMLTag tagServer(*this, "mpi-single", serverOcc, TAG_SERVER);
    doc = "When a solver runs in parallel, it has to use preCICE in form of a ";
    doc += "separatly running server. This is enabled by this tag. ";
    doc += "The communication between participant and server is done by mpi ";
    doc += "with startup in a common communication space.";
    tagServer.setDocumentation(doc);
    serverTags.push_back(tagServer);
  }
  foreach (XMLTag& tagServer, serverTags){
    tag.addSubtag(tagServer);
  }

  std::list<XMLTag> masterTags;
  XMLTag::Occurrence masterOcc = XMLTag::OCCUR_NOT_OR_ONCE;
  {
    XMLTag tagMaster(*this, "sockets", masterOcc, TAG_MASTER);
    doc = "A solver in parallel has to use either a Master or a Server, but not both. ";
    doc += "If you use a Master, you do not have to start-up a further executable, ";
    doc += "all communication is handled peer to peer. One solver process becomes the ";
    doc += " Master handling the synchronization of all slaves. Here, you define then ";
    doc += " the communication between the Master and all slaves. ";
    doc += "The communication between Master and slaves is done by sockets.";
    tagMaster.setDocumentation(doc);

    XMLAttribute<int> attrPort("port");
    doc = "Port number (16-bit unsigned integer) to be used for socket ";
    doc += "communiation. The default is \"0\", what means that OS will ";
    doc += "dynamically search for a free port (if at least one exists) and ";
    doc += "bind it automatically.";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(0);
    tagMaster.addAttribute(attrPort);

    XMLAttribute<std::string> attrNetwork(ATTR_NETWORK);
    doc = "Network name to be used for socket communiation. ";
    doc += "Default is \"lo\", i.e., the local host loopback.";
    attrNetwork.setDocumentation(doc);
    attrNetwork.setDefaultValue("lo");
    tagMaster.addAttribute(attrNetwork);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where connection information is exchanged. By default, the ";
    doc += "directory of startup is chosen.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tagMaster.addAttribute(attrExchangeDirectory);

    masterTags.push_back(tagMaster);
  }
  {
    XMLTag tagMaster(*this, "mpi", masterOcc, TAG_MASTER);
    doc = "A solver in parallel has to use either a Master or a Server, but not both. ";
    doc += "If you use a Master, you do not have to start-up a further executable, ";
    doc += "all communication is handled peer to peer. One solver process becomes the ";
    doc += " Master handling the synchronization of all slaves. Here, you define then ";
    doc += " the communication between the Master and all slaves. ";
    doc += "The communication between Master and slaves is done by mpi ";
    doc += "with startup in separated communication spaces.";
    tagMaster.setDocumentation(doc);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where connection information is exchanged. By default, the ";
    doc += "directory of startup is chosen.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tagMaster.addAttribute(attrExchangeDirectory);

    masterTags.push_back(tagMaster);
  }
  {
    XMLTag tagMaster(*this, "mpi-single", masterOcc, TAG_MASTER);
    doc = "A solver in parallel has to use either a Master or a Server, but not both. ";
    doc += "If you use a Master, you do not have to start-up a further executable, ";
    doc += "all communication is handled peer to peer. One solver process becomes the ";
    doc += " Master handling the synchronization of all slaves. Here, you define then ";
    doc += " the communication between the Master and all slaves. ";
    doc += "The communication between Master and slaves is done by mpi ";
    doc += "with startup in one communication spaces.";
    tagMaster.setDocumentation(doc);

    masterTags.push_back(tagMaster);
  }
  foreach (XMLTag& tagMaster, masterTags){
    tag.addSubtag(tagMaster);
  }

//  XMLAttribute<std::string> attrCom(ATTR_COMMUNICATION);
//  doc = "Sets the communication means to transmit data between solver processes ";
//  doc += "and preCICE server (see tag <communication>).";
//  ValidatorEquals<std::string> validMPI("mpi");
//  ValidatorEquals<std::string> validMPISingle("mpi-single");
//  ValidatorEquals<std::string> validSockets("sockets");
//  attrCom.setValidator(validMPI || validMPISingle || validSockets);
//  attrCom.setDocumentation(doc);
//  tagServer.addAttribute(attrCom);
//  XMLAttribute<std::string> attrContext(ATTR_CONTEXT);
//  doc = "Sets the communication context (see tag <communication>).";
//  attrContext.setDocumentation(doc);
//  attrContext.setDefaultValue("");
//  tagServer.addAttribute(attrContext);
//  tag.addSubtag(tagServer);

  parent.addSubtag(tag);
}

void ParticipantConfiguration:: setDimensions
(
  int dimensions )
{
  preciceTrace1("setDimensions()", dimensions);
  assertion1((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

//bool ParticipantConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  using utils::XMLTag;
//  using utils::XMLAttribute;
//  using utils::ValidatorEquals;
//  XMLTag tagParticipant ( TAG, XMLTag::OCCUR_ONCE );
//
//  XMLAttribute<std::string> attrName ( ATTR_NAME );
//  tagParticipant.addAttribute ( attrName );
//
//  XMLTag tagWriteData ( TAG_WRITE, XMLTag::OCCUR_ARBITRARY);
//  XMLTag tagReadData  ( TAG_READ, XMLTag::OCCUR_ARBITRARY);
//  XMLAttribute<std::string> attrDataName ( ATTR_NAME );
//  tagWriteData.addAttribute ( attrDataName );
//  tagReadData.addAttribute ( attrDataName );
//  XMLAttribute<std::string> attrMesh ( ATTR_MESH );
//  tagWriteData.addAttribute ( attrMesh );
//  tagReadData.addAttribute ( attrMesh );
//  tagParticipant.addSubtag ( tagWriteData );
//  tagParticipant.addSubtag ( tagReadData );
//
//  XMLTag tagMapping ( mapping::MappingConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
//  tagParticipant.addSubtag ( tagMapping );
//
//  XMLTag tagAction ( action::ActionConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
//  tagParticipant.addSubtag ( tagAction );
//
//  XMLTag tagExport ( io::ExportConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
//  tagParticipant.addSubtag ( tagExport );
//
//  XMLTag tagWatchPoint ( TAG_WATCH_POINT, XMLTag::OCCUR_ARBITRARY );
//  tagWatchPoint.addAttribute ( attrName );
//  XMLAttribute<std::string> attrGeometry ( ATTR_MESH );
//  tagWatchPoint.addAttribute ( attrGeometry );
//  if (_dimensions == 2){
//    XMLAttribute<utils::Vector2D> attrCoordinate ( ATTR_COORDINATE );
//    tagWatchPoint.addAttribute ( attrCoordinate );
//  }
//  else {
//    XMLAttribute<utils::Vector3D> attrCoordinate ( ATTR_COORDINATE );
//    tagWatchPoint.addAttribute ( attrCoordinate );
//  }
//  tagParticipant.addSubtag ( tagWatchPoint );
//
//  XMLTag tagUseMesh ( TAG_USE_MESH, XMLTag::OCCUR_ARBITRARY );
//  tagUseMesh.addAttribute ( attrName );
//  if ( _dimensions == 2 ){
//    XMLAttribute<utils::Vector2D> attrLocalOffset ( ATTR_LOCAL_OFFSET );
//    attrLocalOffset.setDefaultValue ( utils::Vector2D(0.0) );
//    tagUseMesh.addAttribute ( attrLocalOffset );
//  }
//  else {
//    XMLAttribute<utils::Vector3D> attrLocalOffset ( ATTR_LOCAL_OFFSET );
//    attrLocalOffset.setDefaultValue ( utils::Vector3D(0.0) );
//    tagUseMesh.addAttribute ( attrLocalOffset );
//  }
//  XMLAttribute<std::string> attrFrom ( ATTR_FROM );
//  attrFrom.setDefaultValue ( "" );
//  tagUseMesh.addAttribute ( attrFrom );
//  XMLAttribute<bool> attrProvide ( ATTR_PROVIDE );
//  attrProvide.setDefaultValue ( false );
//  tagUseMesh.addAttribute ( attrProvide );
//  tagParticipant.addSubtag ( tagUseMesh );
//
//  XMLTag tagServer ( TAG_SERVER, XMLTag::OCCUR_NOT_OR_ONCE );
//  XMLAttribute<std::string> attrCom ( ATTR_COMMUNICATION );
//  tagServer.addAttribute ( attrCom );
//  XMLAttribute<std::string> attrContext(ATTR_CONTEXT);
//  attrContext.setDefaultValue("");
//  tagServer.addAttribute(attrContext);
//  tagParticipant.addSubtag ( tagServer );
//
//  _isValid = tagParticipant.parse ( xmlReader, *this );
//  return _isValid;
//}

void ParticipantConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1 ( "xmlTagCallback()", tag.getName() );
  if (tag.getName() == TAG){
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    impl::PtrParticipant p(new impl::Participant(name, _meshConfig));
    _participants.push_back(p);
  }
  else if (tag.getName() == TAG_USE_MESH){
    assertion(_dimensions != 0); // setDimensions() has been called
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    utils::DynVector offset(_dimensions);
    offset = tag.getDynVectorAttributeValue(ATTR_LOCAL_OFFSET, _dimensions);
    std::string from = tag.getStringAttributeValue(ATTR_FROM);
    double safetyFactor = tag.getDoubleAttributeValue(ATTR_SAFETY_FACTOR);
    if (safetyFactor < 0){
      std::ostringstream stream;
      stream << "Safety Factor must be positive or 0";
      throw stream.str();
    }
    bool provide = tag.getBooleanAttributeValue(ATTR_PROVIDE);
    mesh::PtrMesh mesh = _meshConfig->getMesh(name);
    if (mesh.get() == 0){
      std::ostringstream stream;
      stream << "Participant \"" << _participants.back()->getName()
             << "\" uses mesh \"" << name << "\" which is not defined";
      throw stream.str();
    }
    geometry::PtrGeometry geo ( _geometryConfig->getGeometry(name) );
    spacetree::PtrSpacetree spacetree;
    if (_meshConfig->doesMeshUseSpacetree(name)){
      std::string spacetreeName = _meshConfig->getSpacetreeName(name);
      spacetree = _spacetreeConfig->getSpacetree ( spacetreeName );
    }
    _participants.back()->useMesh ( mesh, geo, spacetree, offset, false, from, safetyFactor, provide );
  }
//  else if ( tag.getName() == mapping::MappingConfiguration::TAG ) {
//    return _mappingConfig->parseSubtag ( xmlReader );
//  }
  else if ( tag.getName() == TAG_WRITE ) {
    std::string dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh mesh = _meshConfig->getMesh ( meshName );
    preciceCheck ( mesh.use_count() > 0, "xmlTagCallback()", "Participant "
                   << "\"" << _participants.back()->getName() << "\" has to use "
                   << "mesh \"" << meshName << "\" in order to write data to it!" );
    mesh::PtrData data = getData ( mesh, dataName );
    _participants.back()->addWriteData ( data, mesh );
  }
  else if ( tag.getName() == TAG_READ ) {
    std::string dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh mesh = _meshConfig->getMesh ( meshName );
    preciceCheck ( mesh.use_count() > 0, "xmlTagCallback()", "Participant "
                   << "\"" << _participants.back()->getName() << "\" has to use "
                   << "mesh \"" << meshName << "\" in order to read data from it!" );
    mesh::PtrData data = getData ( mesh, dataName );
    _participants.back()->addReadData ( data, mesh );
  }
//  else if ( tag.getName() == action::ActionConfiguration::TAG ){
//    action::ActionConfiguration config ( _meshConfig );
//    bool success = config.parseSubtag ( xmlReader );
//    if ( success ){
//      bool used = _participants.back()->isMeshUsed( config.getUsedMeshID() );
//      preciceCheck ( used, "xmlTagCallback()", "Data action of participant "
//                     << _participants.back()->getName()
//                     << "\" uses mesh which is not used by the participant!" );
//      _participants.back()->addAction ( config.getAction() );
//    }
//    return success;
//  }
//  else if ( tag.getName() == io::ExportConfiguration::TAG ){
//    io::ExportConfiguration config;
//    if ( config.parseSubtag(xmlReader) ){
//      //preciceDebug ( "Setting export context with "
//      //               << config.exports().size() << " exports" );
//      _participants.back()->addExportContext ( config.getExportContext() );
//      return true;
//    }
//    return false;
//  }
  else if ( tag.getName() == TAG_WATCH_POINT ){
    assertion(_dimensions != 0); // setDimensions() has been called
    WatchPointConfig config;
    config.name = tag.getStringAttributeValue(ATTR_NAME);
    config.nameMesh = tag.getStringAttributeValue(ATTR_MESH);
    config.coordinates.append(
        tag.getDynVectorAttributeValue(ATTR_COORDINATE, _dimensions));
    _watchPointConfigs.push_back(config);
  }
  else if (tag.getNamespace() == TAG_SERVER){
    com::CommunicationConfiguration comConfig;
    com::Communication::SharedPointer com = comConfig.createCommunication(tag);
    _participants.back()->setClientServerCommunication(com);
  }
  else if (tag.getNamespace() == TAG_MASTER){
    com::CommunicationConfiguration comConfig;
    com::Communication::SharedPointer com = comConfig.createCommunication(tag);
    _participants.back()->setMasterSlaveCommunication(com);
  }
}

void ParticipantConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  if (tag.getName() == TAG){
    finishParticipantConfiguration(_participants.back());
  }
//  else if (tag.getName() == action::ActionConfiguration::TAG){
//    bool used = _participants.back()->isMeshUsed(_actionConfig->getUsedMeshID());
//    preciceCheck(used, "xmlTagCallback()", "Data action of participant "
//                 << _participants.back()->getName()
//                 << "\" uses mesh which is not used by the participant!");
//    _participants.back()->addAction(_actionConfig->getAction());
//  }
}


//bool ParticipantConfiguration:: isValid() const
//{
//  return _isValid;
//}

void ParticipantConfiguration:: addParticipant
(
  const impl::PtrParticipant&             participant,
  const mapping::PtrMappingConfiguration& mappingConfig )
{
  _participants.push_back ( participant );
  _mappingConfig = mappingConfig;
  finishParticipantConfiguration ( participant );
}

const std::vector<impl::PtrParticipant>&
ParticipantConfiguration:: getParticipants() const
{
  return _participants;
}

mesh::PtrMesh ParticipantConfiguration:: copy
(
  const mesh::PtrMesh& mesh ) const
{
  int dim = mesh->getDimensions();
  std::string name(mesh->getName());
  bool flipNormals = mesh->isFlipNormals();
  mesh::Mesh* meshCopy = new mesh::Mesh("Local_" + name, dim, flipNormals);
  foreach (const mesh::PtrData& data, mesh->data()){
    meshCopy->createData(data->getName(), data->getDimensions());
  }
  return mesh::PtrMesh(meshCopy);
}

const mesh::PtrData & ParticipantConfiguration:: getData
(
  const mesh::PtrMesh& mesh,
  const std::string&   nameData ) const
{
  foreach ( const mesh::PtrData & data, mesh->data() ){
    if ( data->getName() == nameData ) {
      return data;
    }
  }
  preciceError ( "getData()", "Participant \"" << _participants.back()->getName()
                 << "\" assignes data \"" << nameData << "\" wrongly to mesh \""
                 << mesh->getName() << "\"!" );
}

void ParticipantConfiguration:: finishParticipantConfiguration
(
  const impl::PtrParticipant& participant )
{
  preciceTrace1("finishParticipantConfiguration()", participant->getName());

  // Set input/output meshes for data mappings and mesh requirements
  typedef mapping::MappingConfiguration::ConfiguredMapping ConfMapping;
  foreach (const ConfMapping& confMapping, _mappingConfig->mappings()){
    int fromMeshID = confMapping.fromMesh->getID();
    int toMeshID = confMapping.toMesh->getID();

    preciceCheck(participant->isMeshUsed(fromMeshID), "finishParticipantConfiguration()",
        "Participant \"" << participant->getName() << "\" has mapping"
        << " from mesh \"" << confMapping.fromMesh->getName() << "\" which he does not use!");
    preciceCheck(participant->isMeshUsed(toMeshID), "finishParticipantConfiguration()",
            "Participant \"" << participant->getName() << "\" has mapping"
            << " to mesh \"" << confMapping.toMesh->getName() << "\" which he does not use!");
    if(participant->useMaster()){
      if((confMapping.direction == mapping::MappingConfiguration::WRITE &&
          confMapping.mapping->getConstraint()==mapping::Mapping::CONSISTENT) ||
         (confMapping.direction == mapping::MappingConfiguration::READ &&
          confMapping.mapping->getConstraint()==mapping::Mapping::CONSERVATIVE)){
        preciceError ( "finishParticipantConfiguration()",
                       "If a participant uses a master parallelization, only the mapping"
                    << " combinations read-consistent and write-conservative are allowed");
      }
    }



    impl::MeshContext& fromMeshContext = participant->meshContext(fromMeshID);
    impl::MeshContext& toMeshContext = participant->meshContext(toMeshID);

    precice::impl::MappingContext* mappingContext = new precice::impl::MappingContext();
    mappingContext->fromMeshID = fromMeshID;
    mappingContext->toMeshID = toMeshID;
    mappingContext->timing = confMapping.timing;

    mapping::PtrMapping& map = mappingContext->mapping;
    assertion(map.get() == NULL);
    map = confMapping.mapping;

    const mesh::PtrMesh& input = fromMeshContext.mesh;
    const mesh::PtrMesh& output = toMeshContext.mesh;
    preciceDebug("Configure mapping for input=" << input->getName()
           << ", output=" << output->getName());
    map->setMeshes(input, output);

    if (confMapping.direction == mapping::MappingConfiguration::WRITE){
      participant->addWriteMappingContext(mappingContext);
    }
    else {
      assertion(confMapping.direction == mapping::MappingConfiguration::READ);
      participant->addReadMappingContext(mappingContext);
    }

    if (map->getInputRequirement() > fromMeshContext.meshRequirement){
      fromMeshContext.meshRequirement = map->getInputRequirement();
    }
    if (map->getOutputRequirement() > toMeshContext.meshRequirement){
      toMeshContext.meshRequirement = map->getOutputRequirement();
    }

    fromMeshContext.fromMappingContext = *mappingContext;
    toMeshContext.toMappingContext = *mappingContext;
  }
  _mappingConfig->resetMappings();

  // Set participant data for data contexts
  foreach (impl::DataContext& dataContext, participant->writeDataContexts()){
    int fromMeshID = dataContext.mesh->getID();
    preciceCheck(participant->isMeshUsed(fromMeshID), "finishParticipant()",
        "Participant \"" << participant->getName() << "\" has to use mesh \""
        << dataContext.mesh->getName() << "\" when writing data to it!");

    foreach (impl::MappingContext& mappingContext, participant->writeMappingContexts()){
      if(mappingContext.fromMeshID==fromMeshID){
        dataContext.mappingContext = mappingContext;
        impl::MeshContext& meshContext = participant->meshContext(mappingContext.toMeshID);
        foreach (mesh::PtrData data, meshContext.mesh->data()){
          if(data->getName()==dataContext.fromData->getName()){
            dataContext.toData = data;
          }
        }
        preciceCheck(dataContext.fromData!=dataContext.toData,"finishParticipant()",
              "The mesh \"" << meshContext.mesh->getName() << "\" needs to use the data \""
              << dataContext.fromData->getName() << "\"! to allow the write mapping");
      }
    }
  }

  foreach (impl::DataContext& dataContext, participant->readDataContexts()){
    int toMeshID = dataContext.mesh->getID();
    preciceCheck(participant->isMeshUsed(toMeshID), "finishParticipant()",
      "Participant \"" << participant->getName() << "\" has to use mesh \""
      << dataContext.mesh->getName() << "\" when writing data to it!");

    foreach (impl::MappingContext& mappingContext, participant->readMappingContexts()){
      if(mappingContext.toMeshID==toMeshID){
        dataContext.mappingContext = mappingContext;
        impl::MeshContext& meshContext = participant->meshContext(mappingContext.fromMeshID);
        foreach (mesh::PtrData data, meshContext.mesh->data()){
          if(data->getName()==dataContext.toData->getName()){
            dataContext.fromData = data;
          }
        }
        preciceCheck(dataContext.toData!=dataContext.fromData,"finishParticipant()",
              "The mesh \"" << meshContext.mesh->getName() << "\" needs to use the data \""
              << dataContext.toData->getName() << "\"! to allow the read mapping");
      }
    }
  }

  // Add actions
  foreach (const action::PtrAction& action, _actionConfig->actions()){
    bool used = _participants.back()->isMeshUsed(action->getMesh()->getID());
    preciceCheck(used, "finishParticipantConfiguration()", "Data action of participant "
                 << _participants.back()->getName()
                 << "\" uses mesh which is not used by the participant!");
    _participants.back()->addAction(action);
  }
  _actionConfig->resetActions();

  // Add export contexts
  foreach (const io::ExportContext& context, _exportConfig->exportContexts()){
    preciceCheck(not participant->useMaster(), "finishParticipantConfiguration()",
        "To use exports while using a master is not yet supported");
    _participants.back()->addExportContext(context);
  }
  _exportConfig->resetExports();

  // Create watch points
  foreach ( const WatchPointConfig & config, _watchPointConfigs ){
    mesh::PtrMesh mesh;
    for ( const impl::MeshContext* context : participant->usedMeshContexts() ){
      if ( context->mesh->getName() == config.nameMesh ){
        mesh = context->mesh;
      }
    }
    preciceCheck ( mesh.use_count() > 0, "xmlEndTagCallback()",
                   "Participant \"" << participant->getName()
                   << "\" defines watchpoint \"" << config.name
                   << "\" for mesh \"" << config.nameMesh
                   << "\" which is not used by him!" );
    std::string filename = config.name + ".watchpoint.txt";
    impl::PtrWatchPoint watchPoint (
        new impl::WatchPoint(config.coordinates, mesh, filename) );
    participant->addWatchPoint ( watchPoint );
  }
  _watchPointConfigs.clear ();
}


}} // namespace precice, config

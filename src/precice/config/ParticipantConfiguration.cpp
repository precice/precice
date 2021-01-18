#include "ParticipantConfiguration.hpp"
#include <algorithm>
#include <list>
#include <memory>
#include <stdexcept>
#include "action/Action.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "io/ExportContext.hpp"
#include "io/ExportVTK.hpp"
#include "io/ExportVTKXML.hpp"
#include "io/SharedPointer.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/MappingContext.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/WatchIntegral.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/PointerVector.hpp"
#include "utils/assertion.hpp"
#include "utils/networking.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace config {

ParticipantConfiguration::ParticipantConfiguration(
    xml::XMLTag &                     parent,
    const mesh::PtrMeshConfiguration &meshConfiguration)
    : _meshConfig(meshConfiguration)
{
  PRECICE_ASSERT(_meshConfig);
  using namespace xml;
  std::string doc;
  XMLTag      tag(*this, TAG, XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Represents one solver using preCICE. At least two ";
  doc += "participants have to be defined.";
  tag.setDocumentation(doc);

  auto attrName = XMLAttribute<std::string>(ATTR_NAME)
                      .setDocumentation(
                          "Name of the participant. Has to match the name given on construction "
                          "of the precice::SolverInterface object used by the participant.");
  tag.addAttribute(attrName);

  XMLTag tagWriteData(*this, TAG_WRITE, XMLTag::OCCUR_ARBITRARY);
  doc = "Sets data to be written by the participant to preCICE. ";
  doc += "Data is defined by using the <data> tag.";
  tagWriteData.setDocumentation(doc);
  XMLTag tagReadData(*this, TAG_READ, XMLTag::OCCUR_ARBITRARY);
  doc = "Sets data to be read by the participant from preCICE. ";
  doc += "Data is defined by using the <data> tag.";
  tagReadData.setDocumentation(doc);
  auto attrDataName = XMLAttribute<std::string>(ATTR_NAME)
                          .setDocumentation("Name of the data.");
  tagWriteData.addAttribute(attrDataName);
  tagReadData.addAttribute(attrDataName);
  auto attrMesh = XMLAttribute<std::string>(ATTR_MESH)
                      .setDocumentation(
                          "Mesh the data belongs to. If data should be read/written to several "
                          "meshes, this has to be specified separately for each mesh.");
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
  auto attrCoordinate = XMLAttribute<Eigen::VectorXd>(ATTR_COORDINATE)
                            .setDocumentation(
                                "The coordinates of the watch point. If the watch point is not put exactly "
                                "on the mesh to observe, the closest projection of the point onto the "
                                "mesh is considered instead, and values/coordinates are interpolated "
                                "linearly to that point.");
  tagWatchPoint.addAttribute(attrCoordinate);
  tag.addSubtag(tagWatchPoint);

  auto attrScaleWitConn = XMLAttribute<bool>(ATTR_SCALE_WITH_CONN)
                              .setDocumentation("Whether the vertex data is scaled with the element area before "
                                                "summing up or not. In 2D, vertex data is scaled with the average length of "
                                                "neighboring edges. In 3D, vertex data is scaled with the average surface of "
                                                "neighboring triangles. If false, vertex data is directly summed up.");
  XMLTag tagWatchIntegral(*this, TAG_WATCH_INTEGRAL, XMLTag::OCCUR_ARBITRARY);
  doc = "A watch integral can be used to follow the transient change of integral data ";
  doc += "and surface area for a given coupling mesh.";
  tagWatchIntegral.setDocumentation(doc);
  doc = "Name of the watch integral. Is taken in combination with the participant ";
  doc += "name to construct the filename the watch integral data is written to.";
  attrName.setDocumentation(doc);
  tagWatchIntegral.addAttribute(attrName);
  doc = "Mesh to be watched.";
  attrMesh.setDocumentation(doc);
  tagWatchIntegral.addAttribute(attrMesh);
  tagWatchIntegral.addAttribute(attrScaleWitConn);
  tag.addSubtag(tagWatchIntegral);

  XMLTag tagUseMesh(*this, TAG_USE_MESH, XMLTag::OCCUR_ARBITRARY);
  doc = "Makes a mesh (see tag <mesh> available to a participant.";
  tagUseMesh.setDocumentation(doc);
  attrName.setDocumentation("Name of the mesh.");
  tagUseMesh.addAttribute(attrName);
  //  XMLAttribute<Eigen::VectorXd> attrLocalOffset(ATTR_LOCAL_OFFSET);
  //  doc = "The mesh can have an offset only applied for the local participant. ";
  //  doc += "Vector-valued example: '1.0; 0.0; 0.0'";
  //  attrLocalOffset.setDocumentation(doc);
  //  attrLocalOffset.setDefaultValue(Eigen::VectorXd::Constant(3, 0));
  //  tagUseMesh.addAttribute(attrLocalOffset);

  auto attrFrom = XMLAttribute<std::string>(ATTR_FROM, "")
                      .setDocumentation(
                          "If a created mesh should be used by "
                          "another solver, this attribute has to specify the creating participant's"
                          " name. The creator has to use the attribute \"provide\" to signal he is "
                          "providing the mesh geometry.");
  tagUseMesh.addAttribute(attrFrom);
  auto attrSafetyFactor = makeXMLAttribute(ATTR_SAFETY_FACTOR, 0.5)
                              .setDocumentation(
                                  "If a mesh is received from another partipant (see tag <from>), it needs to be"
                                  "decomposed at the receiving participant. To speed up this process, "
                                  "a geometric filter (see tag <geometric-filter>), i.e. filtering by bounding boxes around the local mesh, can be used. "
                                  "This safety factor defines by which factor this local information is "
                                  "increased. An example: 0.5 means that the bounding box is 150% of its original size.");
  tagUseMesh.addAttribute(attrSafetyFactor);

  auto attrGeoFilter = XMLAttribute<std::string>(ATTR_GEOMETRIC_FILTER)
                           .setDocumentation(
                               "If a mesh is received from another partipant (see tag <from>), it needs to be"
                               "decomposed at the receiving participant. To speed up this process, "
                               "a geometric filter, i.e. filtering by bounding boxes around the local mesh, can be used. "
                               "Two different variants are implemented: a filter \"on-master\" strategy, "
                               "which is beneficial for a huge mesh and a low number of processors, and a filter "
                               "\"on-slaves\" strategy, which performs better for a very high number of "
                               "processors. Both result in the same distribution (if the safety factor is sufficiently large). "
                               "\"on-master\" is not supported if you use two-level initialization. "
                               "For very asymmetric cases, the filter can also be switched off completely (\"no-filter\").")
                           .setOptions({VALUE_FILTER_ON_MASTER, VALUE_FILTER_ON_SLAVES, VALUE_NO_FILTER})
                           .setDefaultValue(VALUE_FILTER_ON_SLAVES);
  tagUseMesh.addAttribute(attrGeoFilter);

  auto attrProvide = makeXMLAttribute(ATTR_PROVIDE, false)
                         .setDocumentation(
                             "If this attribute is set to \"on\", the "
                             "participant has to create the mesh geometry before initializing preCICE.");
  tagUseMesh.addAttribute(attrProvide);
  tag.addSubtag(tagUseMesh);

  std::list<XMLTag>  masterTags;
  XMLTag::Occurrence masterOcc = XMLTag::OCCUR_NOT_OR_ONCE;
  {
    XMLTag tagMaster(*this, "sockets", masterOcc, TAG_MASTER);
    doc = "A solver in parallel needs a communication between its ranks. ";
    doc += "By default, the participant's MPI_COM_WORLD is reused.";
    doc += "Use this tag to use TCP/IP sockets instead.";
    tagMaster.setDocumentation(doc);

    auto attrPort = makeXMLAttribute("port", 0)
                        .setDocumentation(
                            "Port number (16-bit unsigned integer) to be used for socket "
                            "communiation. The default is \"0\", what means that OS will "
                            "dynamically search for a free port (if at least one exists) and "
                            "bind it automatically.");
    tagMaster.addAttribute(attrPort);

    auto attrNetwork = makeXMLAttribute(ATTR_NETWORK, utils::networking::loopbackInterfaceName())
                           .setDocumentation(
                               "Interface name to be used for socket communiation. "
                               "Default is the cannonical name of the loopback interface of your platform. "
                               "Might be different on supercomputing systems, e.g. \"ib0\" "
                               "for the InfiniBand on SuperMUC. ");
    tagMaster.addAttribute(attrNetwork);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, "")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen.");
    tagMaster.addAttribute(attrExchangeDirectory);

    masterTags.push_back(tagMaster);
  }
  {
    XMLTag tagMaster(*this, "mpi", masterOcc, TAG_MASTER);
    doc = "A solver in parallel needs a communication between its ranks. ";
    doc += "By default, the participant's MPI_COM_WORLD is reused.";
    doc += "Use this tag to use MPI with separated communication spaces instead instead.";
    tagMaster.setDocumentation(doc);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, "")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen.");
    tagMaster.addAttribute(attrExchangeDirectory);

    masterTags.push_back(tagMaster);
  }
  {
    XMLTag tagMaster(*this, "mpi-single", masterOcc, TAG_MASTER);
    doc = "A solver in parallel needs a communication between its ranks. ";
    doc += "By default (which is this option), the participant's MPI_COM_WORLD is reused.";
    doc += "This tag is only used to ensure backwards compatibility.";
    tagMaster.setDocumentation(doc);

    masterTags.push_back(tagMaster);
  }
  for (XMLTag &tagMaster : masterTags) {
    tag.addSubtag(tagMaster);
  }

  parent.addSubtag(tag);
}

void ParticipantConfiguration::setDimensions(
    int dimensions)
{
  PRECICE_TRACE(dimensions);
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

void ParticipantConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == TAG) {
    std::string          name = tag.getStringAttributeValue(ATTR_NAME);
    impl::PtrParticipant p(new impl::Participant(name, _meshConfig));
    _participants.push_back(p);
  } else if (tag.getName() == TAG_USE_MESH) {
    PRECICE_ASSERT(_dimensions != 0); // setDimensions() has been called
    std::string     name = tag.getStringAttributeValue(ATTR_NAME);
    Eigen::VectorXd offset(_dimensions);
    /// @todo offset currently not supported
    //offset = tag.getEigenVectorXdAttributeValue(ATTR_LOCAL_OFFSET, _dimensions);
    std::string                                   from         = tag.getStringAttributeValue(ATTR_FROM);
    double                                        safetyFactor = tag.getDoubleAttributeValue(ATTR_SAFETY_FACTOR);
    partition::ReceivedPartition::GeometricFilter geoFilter    = getGeoFilter(tag.getStringAttributeValue(ATTR_GEOMETRIC_FILTER));
    PRECICE_CHECK(safetyFactor >= 0, "Participant \"" << context.name << "\" uses mesh \"" << name << "\" with safety-factor=\"" << safetyFactor << "\". Please use a positive or zero safety-factor instead.")
    bool provide = tag.getBooleanAttributeValue(ATTR_PROVIDE);
    if (_participants.back()->getName() == from) {
      PRECICE_CHECK(provide, "Participant \"" << context.name << "\" cannot use mesh \"" << name << "\" from itself. Use the \"from\"-field to specify which participant has to communicate the mesh to \"" << context.name << "\".");
    }
    mesh::PtrMesh mesh = _meshConfig->getMesh(name);
    PRECICE_CHECK(mesh, "Participant \"" << _participants.back()->getName()
                                         << "\" uses mesh \"" << name << "\" which is not defined. "
                                         << "Please check the use-mesh node with name=\"" << name << "\" or define the mesh.");
    if ((geoFilter != partition::ReceivedPartition::GeometricFilter::ON_SLAVES || safetyFactor != 0.5) && from == "") {
      PRECICE_CHECK(false, "Participant \"" << _participants.back()->getName()
                                            << "\" uses mesh \"" << name << "\", which is not received (no \"from\"), but has a geometric-filter and/or"
                                            << " a safety factor defined. Please extend the use-mesh tag as follows: <use-mesh name=\"" << name << "\" from=\"(other participant)\" />");
    }
    _participants.back()->useMesh(mesh, offset, false, from, safetyFactor, provide, geoFilter);
  } else if (tag.getName() == TAG_WRITE) {
    std::string   dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string   meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh mesh     = _meshConfig->getMesh(meshName);
    PRECICE_CHECK(mesh, "Participant "
                            << "\"" << _participants.back()->getName() << "\" has to use "
                            << "mesh \"" << meshName << "\" in order to write data to it. "
                            << "Please add a use-mesh node with name=\"" << meshName << "\".");
    mesh::PtrData data = getData(mesh, dataName);
    _participants.back()->addWriteData(data, mesh);
  } else if (tag.getName() == TAG_READ) {
    std::string   dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string   meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh mesh     = _meshConfig->getMesh(meshName);
    PRECICE_CHECK(mesh, "Participant "
                            << "\"" << _participants.back()->getName() << "\" has to use "
                            << "mesh \"" << meshName << "\" in order to read data from it. "
                            << "Please add a use-mesh node with name=\"" << meshName << "\".");
    mesh::PtrData data = getData(mesh, dataName);
    _participants.back()->addReadData(data, mesh);
  } else if (tag.getName() == TAG_WATCH_POINT) {
    PRECICE_ASSERT(_dimensions != 0); // setDimensions() has been called
    WatchPointConfig config;
    config.name        = tag.getStringAttributeValue(ATTR_NAME);
    config.nameMesh    = tag.getStringAttributeValue(ATTR_MESH);
    config.coordinates = tag.getEigenVectorXdAttributeValue(ATTR_COORDINATE, _dimensions);
    _watchPointConfigs.push_back(config);
  } else if (tag.getName() == TAG_WATCH_INTEGRAL) {
    PRECICE_ASSERT(_dimensions != 0);
    WatchIntegralConfig config;
    config.name        = tag.getStringAttributeValue(ATTR_NAME);
    config.nameMesh    = tag.getStringAttributeValue(ATTR_MESH);
    config.isScalingOn = tag.getBooleanAttributeValue(ATTR_SCALE_WITH_CONN);
    _watchIntegralConfigs.push_back(config);
  } else if (tag.getNamespace() == TAG_MASTER) {
    com::CommunicationConfiguration comConfig;
    com::PtrCommunication           com = comConfig.createCommunication(tag);
    utils::MasterSlave::_communication  = com;
    _isMasterDefined                    = true;
    _participants.back()->setUseMaster(true);
  }
}

void ParticipantConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  if (tag.getName() == TAG) {
    finishParticipantConfiguration(context, _participants.back());
  }
}

const std::vector<impl::PtrParticipant> &
ParticipantConfiguration::getParticipants() const
{
  return _participants;
}

partition::ReceivedPartition::GeometricFilter ParticipantConfiguration::getGeoFilter(const std::string &geoFilter) const
{
  if (geoFilter == VALUE_FILTER_ON_MASTER) {
    return partition::ReceivedPartition::GeometricFilter::ON_MASTER;
  } else if (geoFilter == VALUE_FILTER_ON_SLAVES) {
    return partition::ReceivedPartition::GeometricFilter::ON_SLAVES;
  } else {
    PRECICE_ASSERT(geoFilter == VALUE_NO_FILTER);
    return partition::ReceivedPartition::GeometricFilter::NO_FILTER;
  }
}

/// @todo remove
mesh::PtrMesh ParticipantConfiguration::copy(
    const mesh::PtrMesh &mesh) const
{
  int         dim = mesh->getDimensions();
  std::string name(mesh->getName());
  bool        flipNormals = mesh->isFlipNormals();
  mesh::Mesh *meshCopy    = new mesh::Mesh("Local_" + name, dim, flipNormals, mesh::Mesh::MESH_ID_UNDEFINED);
  for (const mesh::PtrData &data : mesh->data()) {
    meshCopy->createData(data->getName(), data->getDimensions());
  }
  return mesh::PtrMesh(meshCopy);
}

const mesh::PtrData &ParticipantConfiguration::getData(
    const mesh::PtrMesh &mesh,
    const std::string &  nameData) const
{
  for (const mesh::PtrData &data : mesh->data()) {
    if (data->getName() == nameData) {
      return data;
    }
  }
  PRECICE_ERROR("Participant \"" << _participants.back()->getName()
                                 << "\" asks for data \"" << nameData << "\" from mesh \""
                                 << mesh->getName() << "\", but this mesh does not use such data. "
                                 << "Please add a use-data tag with name=\"" << nameData << "\" to this mesh.");
}

void ParticipantConfiguration::finishParticipantConfiguration(
    const xml::ConfigurationContext &context,
    const impl::PtrParticipant &     participant)
{
  PRECICE_TRACE(participant->getName());

  // Set input/output meshes for data mappings and mesh requirements
  using ConfMapping = mapping::MappingConfiguration::ConfiguredMapping;
  for (const ConfMapping &confMapping : _mappingConfig->mappings()) {

    checkIllDefinedMappings(confMapping, participant);

    int fromMeshID = confMapping.fromMesh->getID();
    int toMeshID   = confMapping.toMesh->getID();

    PRECICE_CHECK(participant->isMeshUsed(fromMeshID),
                  "Participant \"" << participant->getName() << "\" has mapping"
                                   << " from mesh \"" << confMapping.fromMesh->getName() << "\", without using this mesh. "
                                   << "Please add a use-mesh tag with name=\"" << confMapping.fromMesh->getName() << "\"");
    PRECICE_CHECK(participant->isMeshUsed(toMeshID),
                  "Participant \"" << participant->getName() << "\" has mapping"
                                   << " to mesh \"" << confMapping.toMesh->getName() << "\", without using this mesh. "
                                   << "Please add a use-mesh tag with name=\"" << confMapping.toMesh->getName() << "\"");
    PRECICE_CHECK((participant->isMeshProvided(fromMeshID) || participant->isMeshProvided(toMeshID)),
                  "Participant \"" << participant->getName() << "\" has mapping"
                                   << " from mesh \"" << confMapping.fromMesh->getName() << "\", "
                                   << " to mesh \"" << confMapping.toMesh->getName() << "\", but neither are provided. "
                                   << "Please mark the mesh provided by this participant by configuring its use-mesh tag with provided=\"true\".");

    if (context.size > 1) {
      if ((confMapping.direction == mapping::MappingConfiguration::WRITE &&
           confMapping.mapping->getConstraint() == mapping::Mapping::CONSISTENT) ||
          (confMapping.direction == mapping::MappingConfiguration::READ &&
           confMapping.mapping->getConstraint() == mapping::Mapping::CONSERVATIVE)) {
        PRECICE_ERROR(
            "For a parallel participant, only the mapping"
            << " combinations read-consistent and write-conservative are allowed");
      }
    }

    impl::MeshContext &fromMeshContext = participant->meshContext(fromMeshID);
    impl::MeshContext &toMeshContext   = participant->meshContext(toMeshID);

    if (confMapping.direction == mapping::MappingConfiguration::READ) {
      PRECICE_CHECK(toMeshContext.provideMesh, "A read mapping of participant \""
                                                   << participant->getName() << "\" needs to map TO a provided mesh. Mesh \""
                                                   << confMapping.toMesh->getName() << "\" is not provided. "
                                                   << "Please add a provide=\"yes\" attribute to the participant's use-mesh tag.");
      PRECICE_CHECK(not fromMeshContext.receiveMeshFrom.empty(), "A read mapping of participant \""
                                                                     << participant->getName()
                                                                     << "\" needs to map FROM a received mesh. Mesh \""
                                                                     << confMapping.fromMesh->getName() << "\" is not received. "
                                                                     << "Please add a from=\"(participant)\" attribute to the participant's use-mesh tag.");
    } else {
      PRECICE_CHECK(fromMeshContext.provideMesh, "A write mapping of participant \""
                                                     << participant->getName() << "\" needs to map FROM a provided mesh. Mesh \""
                                                     << confMapping.fromMesh->getName() << "\" is not provided. "
                                                     << "Please add a provide=\"yes\" attribute to the participant's use-mesh tag.");
      PRECICE_CHECK(not toMeshContext.receiveMeshFrom.empty(), "A write mapping of participant \""
                                                                   << participant->getName()
                                                                   << "\" needs to map TO a received mesh. Mesh \""
                                                                   << confMapping.toMesh->getName() << "\" is not received. "
                                                                   << "Please add a from=\"(participant)\" attribute to the participant's use-mesh tag.");
    }

    if (confMapping.isRBF) {
      fromMeshContext.geoFilter = partition::ReceivedPartition::GeometricFilter::NO_FILTER;
      toMeshContext.geoFilter   = partition::ReceivedPartition::GeometricFilter::NO_FILTER;
    }

    precice::impl::MappingContext *mappingContext = new precice::impl::MappingContext();
    mappingContext->fromMeshID                    = fromMeshID;
    mappingContext->toMeshID                      = toMeshID;
    mappingContext->timing                        = confMapping.timing;

    mapping::PtrMapping &map = mappingContext->mapping;
    PRECICE_ASSERT(map.get() == nullptr);
    map = confMapping.mapping;

    const mesh::PtrMesh &input  = fromMeshContext.mesh;
    const mesh::PtrMesh &output = toMeshContext.mesh;
    PRECICE_DEBUG("Configure mapping for input=" << input->getName()
                                                 << ", output=" << output->getName());
    map->setMeshes(input, output);

    if (confMapping.direction == mapping::MappingConfiguration::WRITE) {
      participant->addWriteMappingContext(mappingContext);
    } else {
      PRECICE_ASSERT(confMapping.direction == mapping::MappingConfiguration::READ);
      participant->addReadMappingContext(mappingContext);
    }

    fromMeshContext.meshRequirement = std::max(
        fromMeshContext.meshRequirement, map->getInputRequirement());
    toMeshContext.meshRequirement = std::max(
        toMeshContext.meshRequirement, map->getOutputRequirement());

    fromMeshContext.fromMappingContexts.push_back(*mappingContext);
    toMeshContext.toMappingContexts.push_back(*mappingContext);
  }
  _mappingConfig->resetMappings();

  // Set participant data for data contexts
  for (impl::DataContext &dataContext : participant->writeDataContexts()) {
    int fromMeshID = dataContext.mesh->getID();
    PRECICE_CHECK(participant->isMeshProvided(fromMeshID),
                  "Participant \"" << participant->getName() << "\" has to use and provide mesh \""
                                   << dataContext.mesh->getName() << "\" to be able to write data to it. "
                                   << "Please add a use-mesh node with name=\"" << dataContext.mesh->getName() << "\" and provide=\"true\".");

    for (impl::MappingContext &mappingContext : participant->writeMappingContexts()) {
      if (mappingContext.fromMeshID == fromMeshID) {
        dataContext.mappingContext     = mappingContext;
        impl::MeshContext &meshContext = participant->meshContext(mappingContext.toMeshID);
        for (mesh::PtrData data : meshContext.mesh->data()) {
          if (data->getName() == dataContext.fromData->getName()) {
            dataContext.toData = data;
          }
        }
        PRECICE_CHECK(dataContext.fromData != dataContext.toData,
                      "Mesh \"" << meshContext.mesh->getName() << "\" needs to use data \""
                                << dataContext.fromData->getName() << "\" to allow a write mapping to it. "
                                << "Please add a use-data node with name=\"" << dataContext.fromData->getName() << "\" to this mesh.");
      }
    }
  }

  for (impl::DataContext &dataContext : participant->readDataContexts()) {
    int toMeshID = dataContext.mesh->getID();
    PRECICE_CHECK(participant->isMeshProvided(toMeshID),
                  "Participant \"" << participant->getName() << "\" has to use and provide mesh \""
                                   << dataContext.mesh->getName() << "\" in order to read data from it. "
                                   << "Please add a use-mesh node with name=\"" << dataContext.mesh->getName() << "\" and provide=\"true\".");

    for (impl::MappingContext &mappingContext : participant->readMappingContexts()) {
      if (mappingContext.toMeshID == toMeshID) {
        dataContext.mappingContext     = mappingContext;
        impl::MeshContext &meshContext = participant->meshContext(mappingContext.fromMeshID);
        for (mesh::PtrData data : meshContext.mesh->data()) {
          if (data->getName() == dataContext.toData->getName()) {
            dataContext.fromData = data;
          }
        }
        PRECICE_CHECK(dataContext.toData != dataContext.fromData,
                      "Mesh \"" << meshContext.mesh->getName() << "\" needs to use data \""
                                << dataContext.toData->getName() << "\" to allow a read mapping to it. "
                                << "Please add a use-data node with name=\"" << dataContext.toData->getName() << "\" to this mesh.");
      }
    }
  }

  // Add actions
  for (const action::PtrAction &action : _actionConfig->actions()) {
    bool used = _participants.back()->isMeshUsed(action->getMesh()->getID());
    PRECICE_CHECK(used, "Data action of participant \""
                            << _participants.back()->getName()
                            << "\" uses mesh \"" << action->getMesh()->getName()
                            << "\", which is not used by the participant. "
                            << "Please add a use-mesh node with name=\"" << action->getMesh()->getName() << "\".");
    _participants.back()->addAction(action);
  }
  _actionConfig->resetActions();

  // Add export contexts
  for (io::ExportContext &exportContext : _exportConfig->exportContexts()) {
    io::PtrExport exporter;
    if (exportContext.type == VALUE_VTK) {
      if (context.size > 1) {
        exporter = io::PtrExport(new io::ExportVTKXML(exportContext.plotNormals));
      } else {
        exporter = io::PtrExport(new io::ExportVTK(exportContext.plotNormals));
      }
    } else {
      PRECICE_ERROR("Participant " << _participants.back()->getName()
                                   << " defines an <export/> tag of unknown type \"" << exportContext.type << "\".");
    }
    exportContext.exporter = exporter;

    _participants.back()->addExportContext(exportContext);
  }
  _exportConfig->resetExports();

  // Create watch points
  for (const WatchPointConfig &config : _watchPointConfigs) {
    const impl::MeshContext *meshContext = participant->usedMeshContextByName(config.nameMesh);

    PRECICE_CHECK(meshContext && meshContext->mesh,
                  "Participant \"" << participant->getName()
                                   << "\" defines watchpoint \"" << config.name
                                   << "\" for mesh \"" << config.nameMesh
                                   << "\" which is not used by the participant. "
                                   << "Please add a use-mesh node with name=\"" << config.nameMesh << "\".");
    PRECICE_CHECK(meshContext->provideMesh,
                  "Participant \"" << participant->getName()
                                   << "\" defines watchpoint \"" << config.name
                                   << "\" for the received mesh \"" << config.nameMesh << "\", which is not allowed. "
                                   << "Please move the watchpoint definition to the participant providing mesh \"" << config.nameMesh << "\".");

    std::string         filename = "precice-" + participant->getName() + "-watchpoint-" + config.name + ".log";
    impl::PtrWatchPoint watchPoint(new impl::WatchPoint(config.coordinates, meshContext->mesh, filename));
    participant->addWatchPoint(watchPoint);
  }
  _watchPointConfigs.clear();

  // Create watch integrals
  for (const WatchIntegralConfig &config : _watchIntegralConfigs) {
    const impl::MeshContext *meshContext = participant->usedMeshContextByName(config.nameMesh);

    PRECICE_CHECK(meshContext && meshContext->mesh,
                  "Participant \"" << participant->getName()
                                   << "\" defines watch integral \"" << config.name
                                   << "\" for mesh \"" << config.nameMesh
                                   << "\" which is not used by the participant. "
                                   << "Please add a use-mesh node with name=\"" << config.nameMesh << "\".");
    PRECICE_CHECK(meshContext->provideMesh,
                  "Participant \"" << participant->getName()
                                   << "\" defines watch integral \"" << config.name
                                   << "\" for the received mesh \"" << config.nameMesh << "\", which is not allowed. "
                                   << "Please move the watchpoint definition to the participant providing mesh \"" << config.nameMesh << "\".");

    std::string            filename = "precice-" + participant->getName() + "-watchintegral-" + config.name + ".log";
    impl::PtrWatchIntegral watchIntegral(new impl::WatchIntegral(meshContext->mesh, filename, config.isScalingOn));
    participant->addWatchIntegral(watchIntegral);
  }
  _watchIntegralConfigs.clear();

  // create default master communication if needed
  if (context.size > 1 && not _isMasterDefined && participant->getName() == context.name) {
#ifdef PRECICE_NO_MPI
    PRECICE_ERROR("Implicit master communications for parallel participants are only available if preCICE was built with MPI. "
                  "Either explicitly define a master communication for each parallel participant or rebuild preCICE with \"PRECICE_MPICommunication=ON\".");
#else
    com::PtrCommunication com          = std::make_shared<com::MPIDirectCommunication>();
    utils::MasterSlave::_communication = com;
    participant->setUseMaster(true);
#endif
  }
  _isMasterDefined = false; // to not mess up with previous participant
}

void ParticipantConfiguration::checkIllDefinedMappings(
    const mapping::MappingConfiguration::ConfiguredMapping &mapping,
    const impl::PtrParticipant &                            participant)
{
  PRECICE_TRACE();
  using ConfMapping = mapping::MappingConfiguration::ConfiguredMapping;

  for (const ConfMapping &configuredMapping : _mappingConfig->mappings()) {
    bool sameToMesh   = mapping.toMesh->getName() == configuredMapping.toMesh->getName();
    bool sameFromMesh = mapping.fromMesh->getName() == configuredMapping.fromMesh->getName();
    if (sameToMesh && sameFromMesh) {
      // It's really the same mapping, not a duplicated one. Those are already checked for in MappingConfiguration.
      return;
    }

    if (sameToMesh) {
      for (const mesh::PtrData &data : mapping.fromMesh->data()) {
        for (const mesh::PtrData &configuredData : configuredMapping.fromMesh->data()) {
          bool sameFromData = data->getName() == configuredData->getName();

          if (not sameFromData) {
            continue;
          }

          bool sameDirection = false;

          if (mapping.direction == mapping::MappingConfiguration::WRITE) {
            for (const impl::DataContext &dataContext : participant->writeDataContexts()) {
              sameDirection |= data->getName() == dataContext.getName();
            }
          }
          if (mapping.direction == mapping::MappingConfiguration::READ) {
            for (const impl::DataContext &dataContext : participant->readDataContexts()) {
              sameDirection |= data->getName() == dataContext.getName();
            }
          }
          PRECICE_CHECK(!sameDirection, "There cannot be two mappings to mesh \""
                                            << mapping.toMesh->getName() << "\" "
                                            << "if the meshes from which is mapped contain duplicated data fields "
                                            << "that are also actually mapped on this participant. "
                                            << "Here, both from meshes contain data \"" << data->getName() << "\". "
                                            << "The mapping is not well defined. Which data \"" << data->getName() << "\" "
                                            << "should be mapped to mesh \"" << mapping.toMesh->getName() << "\"?");
        }
      }
    }
  }
}

} // namespace config
} // namespace precice

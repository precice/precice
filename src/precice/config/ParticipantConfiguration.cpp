#include "ParticipantConfiguration.hpp"
#include <algorithm>
#include <boost/range/adaptor/transformed.hpp>
#include <list>
#include <memory>
#include <stdexcept>
#include <utility>

#include "action/Action.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "io/ExportCSV.hpp"
#include "io/ExportContext.hpp"
#include "io/ExportVTK.hpp"
#include "io/ExportVTP.hpp"
#include "io/ExportVTU.hpp"
#include "io/SharedPointer.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/MappingContext.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/WatchIntegral.hpp"
#include "precice/impl/WatchPoint.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"
#include "utils/networking.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::config {

ParticipantConfiguration::ParticipantConfiguration(
    xml::XMLTag &              parent,
    mesh::PtrMeshConfiguration meshConfiguration)
    : _meshConfig(std::move(meshConfiguration))
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
                          "of the precice::Participant object used by the participant.");
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

  _mappingConfig = std::make_shared<mapping::MappingConfiguration>(
      tag, _meshConfig);

  _actionConfig = std::make_shared<action::ActionConfiguration>(
      tag, _meshConfig);

  _exportConfig = std::make_shared<io::ExportConfiguration>(tag);

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

  XMLTag tagProvideMesh(*this, TAG_PROVIDE_MESH, XMLTag::OCCUR_ARBITRARY);
  doc = "Provide a mesh (see tag `<mesh>`) to other participants.";
  tagProvideMesh.setDocumentation(doc);
  attrName.setDocumentation("Name of the mesh to provide.");
  tagProvideMesh.addAttribute(attrName);
  tag.addSubtag(tagProvideMesh);

  XMLTag tagReceiveMesh(*this, TAG_RECEIVE_MESH, XMLTag::OCCUR_ARBITRARY);
  doc = "Makes a remote mesh (see tag `<mesh>`) available to this participant.";
  tagReceiveMesh.setDocumentation(doc);
  attrName.setDocumentation("Name of the mesh to receive.");
  tagReceiveMesh.addAttribute(attrName);
  auto attrFrom = XMLAttribute<std::string>(ATTR_FROM)
                      .setDocumentation("The name of the participant to receive the mesh from. "
                                        "This participant needs to provide the mesh using `<provide-mesh />`.");

  auto attrEnableAccess = makeXMLAttribute(ATTR_API_ACCESS, false)
                              .setDocumentation(
                                  "Enables access to the data on this received mesh via the preCICE API functions without having to map it to a provided mesh. "
                                  "This is required for direct access or just-in-time mappings. "
                                  "A received mesh needs to be decomposed in preCICE using a region of interest, which cannot be inferred, if there are no mappings to or from a provided mesh. "
                                  "In such cases the API function `setMeshAccessRegion()` must be used to define the region of interest. "
                                  "See the user documentation for more information.");
  tagReceiveMesh.addAttribute(attrEnableAccess);
  // @todo: remove with the next breaking release
  auto attrDirectAccess = makeXMLAttribute(ATTR_DIRECT_ACCESS, false)
                              .setDocumentation(
                                  "Deprecated: use \"api-access\" instead.");
  tagReceiveMesh.addAttribute(attrDirectAccess);

  auto attrGeoFilter = XMLAttribute<std::string>(ATTR_GEOMETRIC_FILTER)
                           .setDocumentation(
                               "For parallel execution, a received mesh needs to be decomposed. "
                               "A geometric filter based on bounding-boxes around the local mesh can speed up this process. "
                               "This setting controls if and where this filter is applied. "
                               "`on-primary-rank` is beneficial for a huge mesh and a low number of processors, but is incompatible with two-level initialization. "
                               "`on-secondary-ranks` performs better for a very high number of processors. "
                               "Both result in the same distribution if the safety-factor is sufficiently large. "
                               "`no-filter` may be useful for very asymmetric cases and for debugging. "
                               "If a mapping based on RBFs (rbf-pum,global-rbf) is used, the filter has no influence and is always `no-filter`.")
                           .setOptions({VALUE_NO_FILTER, VALUE_FILTER_ON_PRIMARY_RANK, VALUE_FILTER_ON_SECONDARY_RANKS})
                           .setDefaultValue(VALUE_FILTER_ON_SECONDARY_RANKS);
  tagReceiveMesh.addAttribute(attrGeoFilter);

  tagReceiveMesh.addAttribute(attrFrom);
  auto attrSafetyFactor = makeXMLAttribute(ATTR_SAFETY_FACTOR, 0.5)
                              .setDocumentation(
                                  "The safety factor of the geometric filter uniformly scales the rank-local bounding box by the given factor. "
                                  "A safety-factor of `0.5` means that the bounding box is 150% of its original size.");
  tagReceiveMesh.addAttribute(attrSafetyFactor);

  tag.addSubtag(tagReceiveMesh);

  std::list<XMLTag>  intraCommTags;
  XMLTag::Occurrence intraCommOcc = XMLTag::OCCUR_NOT_OR_ONCE;
  {
    XMLTag tagIntraComm(*this, "sockets", intraCommOcc, TAG_INTRA_COMM);
    doc = "A solver in parallel needs a communication between its ranks. ";
    doc += "By default, the participant's MPI_COM_WORLD is reused.";
    doc += "Use this tag to use TCP/IP sockets instead.";
    tagIntraComm.setDocumentation(doc);

    auto attrPort = makeXMLAttribute("port", 0)
                        .setDocumentation(
                            "Port number (16-bit unsigned integer) to be used for socket "
                            "communication. The default is \"0\", what means that OS will "
                            "dynamically search for a free port (if at least one exists) and "
                            "bind it automatically.");
    tagIntraComm.addAttribute(attrPort);

    auto attrNetwork = makeXMLAttribute(ATTR_NETWORK, utils::networking::loopbackInterfaceName())
                           .setDocumentation(
                               "Interface name to be used for socket communication. "
                               "Default is the canonical name of the loopback interface of your platform. "
                               "Might be different on supercomputing systems, e.g. \"ib0\" "
                               "for the InfiniBand on SuperMUC. ");
    tagIntraComm.addAttribute(attrNetwork);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, ".")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen.");
    tagIntraComm.addAttribute(attrExchangeDirectory);

    intraCommTags.push_back(tagIntraComm);
  }
  {
    XMLTag tagIntraComm(*this, "mpi", intraCommOcc, TAG_INTRA_COMM);
    doc = "A solver in parallel needs a communication between its ranks. ";
    doc += "By default, the participant's MPI_COM_WORLD is reused.";
    doc += "Use this tag to use MPI with separated communication spaces instead instead.";
    tagIntraComm.setDocumentation(doc);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, ".")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen.");
    tagIntraComm.addAttribute(attrExchangeDirectory);

    intraCommTags.push_back(tagIntraComm);
  }

  for (XMLTag &tagIntraComm : intraCommTags) {
    tag.addSubtag(tagIntraComm);
  }
  parent.addSubtag(tag);
}

void ParticipantConfiguration::setExperimental(
    bool experimental)
{
  _experimental = experimental;
  _mappingConfig->setExperimental(_experimental);
}

void ParticipantConfiguration::setRemeshing(
    bool allowed)
{
  _remeshing = allowed;
}

void ParticipantConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == TAG) {
    const std::string &  name = tag.getStringAttributeValue(ATTR_NAME);
    impl::PtrParticipant p(new impl::ParticipantState(name, _meshConfig));
    _participants.push_back(p);
  } else if (tag.getName() == TAG_PROVIDE_MESH) {
    std::string name = tag.getStringAttributeValue(ATTR_NAME);

    mesh::PtrMesh mesh = _meshConfig->getMesh(name);
    PRECICE_CHECK(mesh,
                  R"(Participant "{}" attempts to provide an unknown mesh "{}". <mesh name="{}"> needs to be defined first.)",
                  _participants.back()->getName(), name, name);
    _participants.back()->provideMesh(mesh);
  } else if (tag.getName() == TAG_RECEIVE_MESH) {
    std::string                                   name              = tag.getStringAttributeValue(ATTR_NAME);
    std::string                                   from              = tag.getStringAttributeValue(ATTR_FROM);
    double                                        safetyFactor      = tag.getDoubleAttributeValue(ATTR_SAFETY_FACTOR);
    partition::ReceivedPartition::GeometricFilter geoFilter         = getGeoFilter(tag.getStringAttributeValue(ATTR_GEOMETRIC_FILTER));
    const bool                                    allowDirectAccess = tag.getBooleanAttributeValue(ATTR_API_ACCESS) || tag.getBooleanAttributeValue(ATTR_DIRECT_ACCESS);
    PRECICE_WARN_IF(tag.getBooleanAttributeValue(ATTR_DIRECT_ACCESS), "The 'direct-access' flag (<receive-mesh direct-access=\"...\" />) is deprecated and will be removed in preCICE v4. Use 'api-access' instead (<receive-mesh api-access=\"...\" />).");

    // Start with defining the mesh
    mesh::PtrMesh mesh = _meshConfig->getMesh(name);
    PRECICE_CHECK(mesh,
                  R"(Participant "{}" attempts to provide an unknown mesh "{}". <mesh name="{}"> needs to be defined first.)",
                  _participants.back()->getName(), name, name);

    // Then check the attributes
    PRECICE_CHECK(!from.empty(),
                  R"(Participant "{}" receives mesh "{}", but doesn't specify where from. )"
                  "Please add the name of the other participant to the receive-mesh tag: <receive-mesh name=\"{}\" from=\"(other participant)\" ... />",
                  context.name, name, name);

    PRECICE_CHECK(_participants.back()->getName() != from,
                  "Participant \"{}\" cannot receive mesh \"{}\" from itself. "
                  "To provide a mesh, use <provide-mesh name=\"{}\" /> instead.",
                  context.name, name, name);

    PRECICE_CHECK(safetyFactor >= 0,
                  "Participant \"{}\" receives mesh \"{}\" with safety-factor=\"{}\". "
                  "Please use a positive or zero safety-factor instead.",
                  context.name, name, safetyFactor);

    _participants.back()->receiveMesh(mesh, from, safetyFactor, geoFilter, allowDirectAccess);
  } else if (tag.getName() == TAG_WRITE) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string        meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh      mesh     = _meshConfig->getMesh(meshName);
    PRECICE_CHECK(mesh,
                  R"(Participant "{}" attempts to write data "{}" from an unknown mesh "{}". <mesh name="{}"> needs to be defined first.)",
                  _participants.back()->getName(), dataName, meshName, meshName);
    mesh::PtrData data = getData(mesh, dataName);
    _participants.back()->addWriteData(data, mesh);
  } else if (tag.getName() == TAG_READ) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_NAME);
    std::string        meshName = tag.getStringAttributeValue(ATTR_MESH);
    mesh::PtrMesh      mesh     = _meshConfig->getMesh(meshName);
    PRECICE_CHECK(mesh,
                  R"(Participant "{}" attempts to read data "{}" to an unknown mesh "{}". <mesh name="{}"> needs to be defined first.)",
                  _participants.back()->getName(), dataName, meshName, meshName);
    mesh::PtrData data = getData(mesh, dataName);
    _participants.back()->addReadData(data, mesh);
  } else if (tag.getName() == TAG_WATCH_POINT) {
    WatchPointConfig config;
    config.name        = tag.getStringAttributeValue(ATTR_NAME);
    config.nameMesh    = tag.getStringAttributeValue(ATTR_MESH);
    config.coordinates = tag.getEigenVectorXdAttributeValue(ATTR_COORDINATE);
    _watchPointConfigs.push_back(config);
  } else if (tag.getName() == TAG_WATCH_INTEGRAL) {
    WatchIntegralConfig config;
    config.name        = tag.getStringAttributeValue(ATTR_NAME);
    config.nameMesh    = tag.getStringAttributeValue(ATTR_MESH);
    config.isScalingOn = tag.getBooleanAttributeValue(ATTR_SCALE_WITH_CONN);
    _watchIntegralConfigs.push_back(config);
  } else if (tag.getNamespace() == TAG_INTRA_COMM) {
    com::CommunicationConfiguration comConfig;
    utils::IntraComm::getCommunication() = comConfig.createCommunication(tag);
    _isIntraCommDefined                  = true;
    _participants.back()->setUsePrimaryRank(true);
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

const impl::PtrParticipant ParticipantConfiguration::getParticipant(const std::string &participantName) const
{
  auto participant = std::find_if(_participants.begin(), _participants.end(), [&participantName](const auto &p) { return p->getName() == participantName; });
  PRECICE_ASSERT(participant != _participants.end(), "Did not find participant \"{}\"", participantName);

  return *participant;
}

std::set<std::string> ParticipantConfiguration::knownParticipants() const
{
  auto range = _participants | boost::adaptors::transformed([](auto &p) { return p->getName(); });
  return {range.begin(), range.end()};
}

bool ParticipantConfiguration::hasParticipant(std::string_view name) const
{
  return std::any_of(_participants.begin(), _participants.end(), [name](auto &p) { return p->getName() == name; });
}

std::string ParticipantConfiguration::hintFor(std::string_view wrongName) const
{
  PRECICE_ASSERT(!hasParticipant(wrongName));

  const auto names   = knownParticipants();
  const auto matches = utils::computeMatches(wrongName, names);

  // Typo detection
  if (matches.front().distance < 3) {
    return fmt::format("Did you mean: \"{}\"?", matches.front().name);
  }

  return fmt::format("Available participants are: {}.", fmt::join(names, ", "));
}

partition::ReceivedPartition::GeometricFilter ParticipantConfiguration::getGeoFilter(const std::string &geoFilter) const
{
  if (geoFilter == VALUE_FILTER_ON_PRIMARY_RANK) {
    return partition::ReceivedPartition::GeometricFilter::ON_PRIMARY_RANK;
  } else if (geoFilter == VALUE_FILTER_ON_SECONDARY_RANKS) {
    return partition::ReceivedPartition::GeometricFilter::ON_SECONDARY_RANKS;
  } else {
    PRECICE_ASSERT(geoFilter == VALUE_NO_FILTER);
    return partition::ReceivedPartition::GeometricFilter::NO_FILTER;
  }
}

const mesh::PtrData &ParticipantConfiguration::getData(
    const mesh::PtrMesh &mesh,
    const std::string &  nameData) const
{
  PRECICE_CHECK(mesh->hasDataName(nameData),
                "Participant \"{}\" asks for data \"{}\" from mesh \"{}\", but this mesh does not use such data. "
                "Please add a use-data tag with name=\"{}\" to this mesh.",
                _participants.back()->getName(), nameData, mesh->getName(), nameData);
  return mesh->data(nameData);
}

void ParticipantConfiguration::finishParticipantConfiguration(
    const xml::ConfigurationContext &context,
    const impl::PtrParticipant &     participant)
{
  PRECICE_TRACE(participant->getName());

  // Set input/output meshes for data mappings and mesh requirements
  // This for loop transforms the MappingConfiguration::ConfiguredMappings
  // into a MappingContext
  using ConfMapping = mapping::MappingConfiguration::ConfiguredMapping;
  for (const ConfMapping &confMapping : _mappingConfig->mappings()) {

    checkIllDefinedMappings(confMapping, participant);

    auto fromMesh = confMapping.fromMesh->getName();
    auto toMesh   = confMapping.toMesh->getName();

    // sanity checks
    if (confMapping.direction == mapping::MappingConfiguration::Direction::READ) {
      // A read mapping maps from received to provided
      PRECICE_CHECK(participant->isMeshReceived(fromMesh) || confMapping.toMesh->isJustInTime() || participant->isMeshProvided(toMesh),
                    "A read mapping of participant \"{}\" needs to map from a received to a provided mesh, but in this case they are swapped. "
                    "Did you intent to map from mesh \"{}\" to mesh \"{}\", or use a write mapping instead?",
                    participant->getName(), confMapping.toMesh->getName(), confMapping.fromMesh->getName());
      PRECICE_CHECK(participant->isMeshReceived(fromMesh),
                    "Participant \"{}\" has a read mapping from mesh \"{}\", without receiving it. "
                    "Please add a receive-mesh tag with name=\"{}\"",
                    participant->getName(), fromMesh, fromMesh);
      // The just-in-time mesh cannot be on the "from" mesh, as only the combinations read-consistent and write-conservative are allowed
      PRECICE_CHECK(confMapping.toMesh->isJustInTime() || participant->isMeshProvided(toMesh),
                    "Participant \"{}\" has a read mapping to mesh \"{}\", without providing it. "
                    "Please add a provide-mesh tag with name=\"{}\"",
                    participant->getName(), toMesh, toMesh);
    } else {
      // A write mapping maps from provided to received
      PRECICE_CHECK(confMapping.fromMesh->isJustInTime() || participant->isMeshProvided(fromMesh) || participant->isMeshReceived(toMesh),
                    "A write mapping of participant \"{}\" needs to map from a provided to a received mesh, but in this case they are swapped. "
                    "Did you intent to map from mesh \"{}\" to mesh \"{}\", or use a read mapping instead?",
                    participant->getName(), confMapping.toMesh->getName(), confMapping.fromMesh->getName());
      // The just-in-time mesh cannot be on the "to" mesh, as only the combinations read-consistent and write-conservative are allowed
      PRECICE_CHECK(confMapping.fromMesh->isJustInTime() || participant->isMeshProvided(fromMesh),
                    "Participant \"{}\" has a write mapping from mesh \"{}\", without providing it. "
                    "Please add a provided-mesh tag with name=\"{}\"",
                    participant->getName(), fromMesh, fromMesh);
      PRECICE_CHECK(participant->isMeshReceived(toMesh),
                    "Participant \"{}\" has a write mapping to mesh \"{}\", without receiving it. "
                    "Please add a receive-mesh tag with name=\"{}\"",
                    participant->getName(), toMesh, toMesh);
    }

    if (context.size > 1 && context.name == participant->getName()) {
      if ((confMapping.direction == mapping::MappingConfiguration::WRITE &&
           confMapping.mapping->getConstraint() == mapping::Mapping::CONSISTENT) ||
          (confMapping.direction == mapping::MappingConfiguration::READ &&
           confMapping.mapping->getConstraint() == mapping::Mapping::CONSERVATIVE)) {
        PRECICE_ERROR("For a parallel participant, only the mapping combinations read-consistent and write-conservative are allowed");
      } else if (confMapping.mapping->isScaledConsistent()) {
        PRECICE_ERROR("Scaled consistent mapping is not yet supported for a parallel participant. "
                      "You could run in serial or use a plain (read-)consistent mapping instead.");
      }
    }

    PRECICE_CHECK(!confMapping.mapping->isScaledConsistent() || !(confMapping.fromMesh->isJustInTime() || confMapping.toMesh->isJustInTime()),
                  "The just-in-time mapping from mesh \"{}\" to mesh \"{}\" was configured with a scaled-consistent constraint. A scaled-consistent constraint is not implemented for just-in-time mappings in preCICE.", confMapping.fromMesh->getName(), confMapping.toMesh->getName());

    // We disable the geometric filter for any kernel method, as the default safety factor is not reliable enough to provide a robust
    // safety margin such that the mapping is still correct.
    if (confMapping.requiresBasisFunction) {
      if (!confMapping.fromMesh->isJustInTime()) {
        participant->meshContext(fromMesh).geoFilter = partition::ReceivedPartition::GeometricFilter::NO_FILTER;
      }
      if (!confMapping.toMesh->isJustInTime()) {
        participant->meshContext(toMesh).geoFilter = partition::ReceivedPartition::GeometricFilter::NO_FILTER;
      }
    }

    // Now we create the mappingContext, which will be stored permanently
    precice::impl::MappingContext mappingContext;
    // Copy over data from MappingConfiguration
    // 1. the mesh data
    mappingContext.fromMeshID = confMapping.fromMesh->getID();
    mappingContext.toMeshID   = confMapping.toMesh->getID();

    // Upon creation, the mapping should be empty
    mapping::PtrMapping &map = mappingContext.mapping;
    PRECICE_ASSERT(map.get() == nullptr);
    // 2. ... and the mappings
    map                                   = confMapping.mapping;
    mappingContext.configuredWithAliasTag = confMapping.configuredWithAliasTag;

    // Set input and output meshes in the Mapping from the mesh contexts
    const mesh::PtrMesh &input  = confMapping.fromMesh->isJustInTime() ? confMapping.fromMesh : participant->meshContext(fromMesh).mesh;
    const mesh::PtrMesh &output = confMapping.toMesh->isJustInTime() ? confMapping.toMesh : participant->meshContext(toMesh).mesh;
    PRECICE_DEBUG("Configure mapping for input={}, output={}", input->getName(), output->getName());
    map->setMeshes(input, output);

    // just-in-time mappings go for now into the participant's mapping context
    // Add the mapping context to the participant, separated by direction
    if (confMapping.direction == mapping::MappingConfiguration::WRITE) {
      participant->addWriteMappingContext(mappingContext);
    } else {
      PRECICE_ASSERT(confMapping.direction == mapping::MappingConfiguration::READ);
      participant->addReadMappingContext(mappingContext);
    }

    // configure the involved mesh context with connectivity requirements stemming from the mapping
    // Add the mapping context to the mesh context, only required to later on forward them to the Partition
    if (!input->isJustInTime()) {
      participant->configureInputMeshContext(fromMesh, mappingContext, map->getInputRequirement());
    }
    if (!output->isJustInTime()) {
      participant->configureOutputMeshContext(toMesh, mappingContext, map->getOutputRequirement());
    }
  }
  // clear the data structure we just transformed and don't need anymore
  _mappingConfig->resetMappings();

  // Now we have the MappingContexts and need to add information on the associated Data we want to map
  //
  // First in write direction:
  // for all writeMappingContexts ...
  for (impl::MappingContext &mappingContext : participant->writeMappingContexts()) {
    // Check, whether we can find a corresponding write data context
    bool dataFound = false;
    for (auto &dataContext : participant->writeDataContexts()) {
      // First we look for the "from" mesh ID from the "data perspective"
      const int fromMeshID = dataContext.getMeshID();
      // and compare it against the "from" mesh ID from the "mapping perspective"
      if (mappingContext.fromMeshID == fromMeshID) {
        // If these two are the same, we have a match of data and mapping contexts on the 'from' side
        //
        // the data context carries now the information about the associated name of the data itself
        // Hence, we look if the "to" mesh (ID) stored in the mappingContext exists on the participant...
        impl::MeshContext &meshContext = participant->meshContext(mappingContext.mapping->getOutputMesh()->getName());
        // .. and if the mesh 'uses' the data to be mapped
        // If this is true, we actually found a proper configuration and add the mapping
        // If it is false, we look for another "from" mesh ID in the data context, because we might have multiple read and write mappings from the same 'from' mesh
        if (meshContext.mesh->hasDataName(dataContext.getDataName())) {
          // Check, if the fromMesh is a provided mesh
          PRECICE_CHECK(participant->isMeshProvided(dataContext.getMeshName()),
                        "Participant \"{}\" has to provide mesh \"{}\" to be able to write data to it. "
                        "Please add a provide-mesh node with name=\"{}\".",
                        participant->getName(), dataContext.getMeshName(), dataContext.getMeshName());
          // here, the mappingContext receives its to and from data pointer
          // we append the mappingContext into the dataContext by copying it over, which is fine, since the context
          // structures operate only on shared object pointers
          dataContext.appendMappingConfiguration(mappingContext, meshContext);
          // Enable gradient data if required
          if (mappingContext.mapping->requiresGradientData() == true) {
            mappingContext.requireGradientData(dataContext.getDataName());
          }
          dataFound = true;
        }
      } else if (mappingContext.mapping->getInputMesh()->isJustInTime()) {
        const int toMeshID = dataContext.getMeshID();
        // We compare here the to mesh instead of the from mesh
        if (mappingContext.toMeshID == toMeshID) {
          impl::MeshContext &meshContext = participant->meshContext(mappingContext.mapping->getOutputMesh()->getName());
          dataContext.addJustInTimeMapping(mappingContext, meshContext);
          if (mappingContext.mapping->requiresGradientData() == true) {
            mappingContext.requireGradientData(dataContext.getDataName());
          }
        }
        dataFound = true;
      }
    }
    PRECICE_CHECK(dataFound,
                  "Participant \"{}\" defines a write mapping from mesh \"{}\" to mesh \"{}\", "
                  "but there is either no corresponding write-data tag or the meshes used "
                  "by this participant lack the necessary use-data tags.",
                  participant->getName(), mappingContext.mapping->getInputMesh()->getName(), mappingContext.mapping->getOutputMesh()->getName());
  }

  // Iterate over all read mappings
  for (impl::MappingContext &mappingContext : participant->readMappingContexts()) {
    // Check, weather we can find a corresponding read data context
    bool dataFound = false;
    for (auto &dataContext : participant->readDataContexts()) {
      // First we look for the "to" mesh ID
      const int toMeshID = dataContext.getMeshID();
      if (mappingContext.toMeshID == toMeshID) {
        // Second we look for the "from" mesh ID
        impl::MeshContext &meshContext = participant->meshContext(mappingContext.mapping->getInputMesh()->getName());
        // If this is true, we actually found a proper configuration
        // If it is false, we look for another "from" mesh ID, because we might have multiple read and write mappings
        if (meshContext.mesh->hasDataName(dataContext.getDataName())) {
          // Check, if the toMesh is a provided mesh
          PRECICE_CHECK(participant->isMeshProvided(dataContext.getMeshName()),
                        "Participant \"{}\" has to provide mesh \"{}\" in order to read data from it. "
                        "Please add a provide-mesh node with name=\"{}\".",
                        participant->getName(), dataContext.getMeshName(), dataContext.getMeshName());
          dataContext.appendMappingConfiguration(mappingContext, meshContext);
          // Enable gradient data if required
          if (mappingContext.mapping->requiresGradientData() == true) {
            mappingContext.requireGradientData(dataContext.getDataName());
          }
          dataFound = true;
        }
      } else if (mappingContext.mapping->getOutputMesh()->isJustInTime()) {
        const int fromMeshID = dataContext.getMeshID();
        // We compare here the from mesh instead of the to mesh
        if (mappingContext.fromMeshID == fromMeshID) {
          impl::MeshContext &meshContext = participant->meshContext(mappingContext.mapping->getInputMesh()->getName());

          dataContext.addJustInTimeMapping(mappingContext, meshContext);
          if (mappingContext.mapping->requiresGradientData() == true) {
            mappingContext.requireGradientData(dataContext.getDataName());
          }
        }
        dataFound = true;
      }
    }
    PRECICE_CHECK(dataFound,
                  "Participant \"{}\" defines a read mapping from mesh \"{}\" to mesh \"{}\", "
                  "but there is either no corresponding read-data tag or the meshes used "
                  "by this participant lack the necessary use-data tags.",
                  participant->getName(), mappingContext.mapping->getInputMesh()->getName(), mappingContext.mapping->getOutputMesh()->getName());
  }

  // Add actions
  for (const action::PtrAction &action : _actionConfig->actions()) {
    bool used = _participants.back()->isMeshUsed(action->getMesh()->getName());
    PRECICE_CHECK(used,
                  "Data action of participant \"{}\" uses mesh \"{}\", which is not used by the participant. "
                  "Please add a provide-mesh or receive-mesh node with name=\"{}\".",
                  _participants.back()->getName(), action->getMesh()->getName(), action->getMesh()->getName());
  }
  for (action::PtrAction &action : _actionConfig->extractActions()) {
    _participants.back()->addAction(std::move(action));
  }

  // Check for unsupported remeshing options
  for (auto &context : participant->writeDataContexts()) {
    PRECICE_CHECK(participant->meshContext(context.getMeshName()).provideMesh || !(participant->isDirectAccessAllowed(context.getMeshName()) && _remeshing), "Writing data via API access (configuration <write-data ... mesh=\"{}\") is not (yet) supported with remeshing", context.getMeshName());
  }

  // Add export contexts
  for (io::ExportContext &exportContext : _exportConfig->exportContexts()) {
    auto kind = exportContext.everyIteration ? io::Export::ExportKind::Iterations : io::Export::ExportKind::TimeWindows;
    // Create one exporter per mesh
    for (const auto &meshContext : participant->usedMeshContexts()) {

      exportContext.meshName = meshContext->mesh->getName();

      io::PtrExport exporter;
      if (exportContext.type == VALUE_VTK) {
        // This is handled with respect to the current configuration context.
        // Hence, this is potentially wrong for every participant other than context.name.
        if (context.size > 1) {
          // Only display the warning message if this participant configuration is the current one.
          if (context.name == participant->getName()) {
            PRECICE_ERROR("You attempted to use the legacy VTK exporter with the parallel participant {}, which isn't supported."
                          "Migrate to another exporter, such as the VTU exporter by specifying \"<export:vtu ... />\"  instead of \"<export:vtk ... />\".",
                          participant->getName());
          }
        } else {
          exporter = io::PtrExport(new io::ExportVTK(
              participant->getName(),
              exportContext.location,
              *meshContext->mesh,
              kind,
              exportContext.everyNTimeWindows,
              context.rank,
              context.size));
        }
      } else if (exportContext.type == VALUE_VTU) {
        exporter = io::PtrExport(new io::ExportVTU(
            participant->getName(),
            exportContext.location,
            *meshContext->mesh,
            kind,
            exportContext.everyNTimeWindows,
            context.rank,
            context.size));
      } else if (exportContext.type == VALUE_VTP) {
        exporter = io::PtrExport(new io::ExportVTP(
            participant->getName(),
            exportContext.location,
            *meshContext->mesh,
            kind,
            exportContext.everyNTimeWindows,
            context.rank,
            context.size));
      } else if (exportContext.type == VALUE_CSV) {
        exporter = io::PtrExport(new io::ExportCSV(
            participant->getName(),
            exportContext.location,
            *meshContext->mesh,
            kind,
            exportContext.everyNTimeWindows,
            context.rank,
            context.size));
      } else {
        PRECICE_ERROR("Participant {} defines an <export/> tag of unknown type \"{}\".",
                      _participants.back()->getName(), exportContext.type);
      }
      exportContext.exporter = std::move(exporter);

      _participants.back()->addExportContext(exportContext);
    }
    PRECICE_WARN_IF(exportContext.everyNTimeWindows > 1 && exportContext.everyIteration,
                    "Participant {} defines an exporter of type {} which exports every iteration. "
                    "This overrides the every-n-time-window value you provided.",
                    _participants.back()->getName(), exportContext.type);
  }
  _exportConfig->resetExports();

  // Create watch points
  if (context.name == participant->getName()) {
    for (const WatchPointConfig &config : _watchPointConfigs) {
      PRECICE_CHECK(participant->isMeshUsed(config.nameMesh),
                    "Participant \"{}\" defines watchpoint \"{}\" for mesh \"{}\" which is not provided by the participant. "
                    "Please add <provide-mesh name=\"{}\" /> to the participant.",
                    participant->getName(), config.name, config.nameMesh, config.nameMesh);
      const auto &meshContext = participant->usedMeshContext(config.nameMesh);
      PRECICE_CHECK(meshContext.provideMesh,
                    "Participant \"{}\" defines watchpoint \"{}\" for the received mesh \"{}\", which is not allowed. "
                    "Please move the watchpoint definition to the participant providing mesh \"{}\".",
                    participant->getName(), config.name, config.nameMesh, config.nameMesh);
      PRECICE_CHECK(config.coordinates.size() == meshContext.mesh->getDimensions(),
                    "Provided coordinate to watch is {}D, which does not match the dimension of the {}D mesh \"{}\".",
                    config.coordinates.size(), meshContext.mesh->getDimensions(), meshContext.mesh->getName());
      std::string filename = "precice-" + participant->getName() + "-watchpoint-" + config.name + ".log";
      participant->addWatchPoint(std::make_shared<impl::WatchPoint>(config.coordinates, meshContext.mesh, std::move(filename)));
    }
  }
  _watchPointConfigs.clear();

  // Create watch integrals
  if (context.name == participant->getName()) {
    for (const WatchIntegralConfig &config : _watchIntegralConfigs) {
      PRECICE_CHECK(participant->isMeshUsed(config.nameMesh),
                    "Participant \"{}\" defines watch integral \"{}\" for mesh \"{}\" which is not used by the participant. "
                    "Please add a provide-mesh node with name=\"{}\".",
                    participant->getName(), config.name, config.nameMesh, config.nameMesh);
      const auto &meshContext = participant->usedMeshContext(config.nameMesh);
      PRECICE_CHECK(meshContext.provideMesh,
                    "Participant \"{}\" defines watch integral \"{}\" for the received mesh \"{}\", which is not allowed. "
                    "Please move the watchpoint definition to the participant providing mesh \"{}\".",
                    participant->getName(), config.name, config.nameMesh, config.nameMesh);

      std::string filename = "precice-" + participant->getName() + "-watchintegral-" + config.name + ".log";
      participant->addWatchIntegral(std::make_shared<impl::WatchIntegral>(meshContext.mesh, std::move(filename), config.isScalingOn));
    }
  }
  _watchIntegralConfigs.clear();

  // create default primary communication if needed
  if (context.size > 1 && not _isIntraCommDefined && participant->getName() == context.name) {
#ifdef PRECICE_NO_MPI
    PRECICE_ERROR("Implicit intra-participant communications for parallel participants are only available if preCICE was built with MPI. "
                  "Either explicitly define an intra-participant communication for each parallel participant or rebuild preCICE with \"PRECICE_MPICommunication=ON\".");
#else
    utils::IntraComm::getCommunication() = std::make_shared<com::MPIDirectCommunication>();
    participant->setUsePrimaryRank(true);
#endif
  }
  _isIntraCommDefined = false; // to not mess up with previous participant
}

void ParticipantConfiguration::checkIllDefinedMappings(
    const mapping::MappingConfiguration::ConfiguredMapping &mapping,
    const impl::PtrParticipant &                            participant)
{
  PRECICE_TRACE();
  using ConfMapping = mapping::MappingConfiguration::ConfiguredMapping;

  for (const ConfMapping &configuredMapping : _mappingConfig->mappings()) {
    bool sameToMesh   = (mapping.toMesh->getName() == configuredMapping.toMesh->getName()) && !mapping.toMesh->isJustInTime();
    bool sameFromMesh = (mapping.fromMesh->getName() == configuredMapping.fromMesh->getName()) && !mapping.fromMesh->isJustInTime();
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
            for (const auto &dataContext : participant->writeDataContexts()) {
              sameDirection |= data->getName() == dataContext.getDataName();
            }
          }
          if (mapping.direction == mapping::MappingConfiguration::READ) {
            for (const auto &dataContext : participant->readDataContexts()) {
              sameDirection |= data->getName() == dataContext.getDataName();
            }
          }
          PRECICE_CHECK(!sameDirection,
                        "There cannot be two mappings to mesh \"{}\" if the meshes from which is mapped contain "
                        "duplicated data fields that are also actually mapped on this participant. "
                        "Here, both from meshes contain data \"{}\". "
                        "The mapping is not well defined. "
                        "Which data \"{}\" should be mapped to mesh \"{}\"?",
                        mapping.toMesh->getName(), data->getName(), data->getName(), mapping.toMesh->getName());
        }
      }
    }
  }
}

} // namespace precice::config

#include "SolverInterfaceConfiguration.hpp"
#include "ParticipantConfiguration.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"

namespace precice {
namespace config {

SolverInterfaceConfiguration:: SolverInterfaceConfiguration(xml::XMLTag& parent )
{
  using namespace xml;
  XMLTag tag(*this, "solver-interface", XMLTag::OCCUR_ONCE);
  tag.setDocumentation("Configuration of simulation relevant features.");
  auto attrDimensions = makeXMLAttribute("dimensions", 2)
                             .setDocumentation("Determines the spatial dimensionality of the configuration")
                             .setOptions({2, 3});
  tag.addAttribute(attrDimensions);

  _dataConfiguration = mesh::PtrDataConfiguration (
      new mesh::DataConfiguration(tag) );
  _meshConfiguration = mesh::PtrMeshConfiguration (
      new mesh::MeshConfiguration(tag, _dataConfiguration) );
  _m2nConfiguration = m2n::M2NConfiguration::SharedPointer (
      new m2n::M2NConfiguration(tag) );
  _participantConfiguration = config::PtrParticipantConfiguration (
    new ParticipantConfiguration(tag, _meshConfiguration) );
  _couplingSchemeConfiguration = cplscheme::PtrCouplingSchemeConfiguration (
    new cplscheme::CouplingSchemeConfiguration(tag, _meshConfiguration,
    _m2nConfiguration) );

  parent.addSubtag(tag);
}

void SolverInterfaceConfiguration:: xmlTagCallback
(
  const xml::ConfigurationContext& context,
  xml::XMLTag& tag )
{
  PRECICE_TRACE();
  if (tag.getName() == "solver-interface"){
    _dimensions = tag.getIntAttributeValue("dimensions");
    _dataConfiguration->setDimensions(_dimensions);
    _meshConfiguration->setDimensions(_dimensions);
    _participantConfiguration->setDimensions(_dimensions);
  }
  else {
    PRECICE_ERROR("Received callback from tag " << tag.getName());
  }
}

void SolverInterfaceConfiguration:: xmlEndTagCallback
(
  const xml::ConfigurationContext& context,
  xml::XMLTag& tag )
{
  PRECICE_TRACE();
  if (tag.getName() == "solver-interface"){
    _meshConfiguration->setMeshSubIDs();

    //test if both participants do have the exchange meshes
    typedef std::map<std::string, std::vector<std::string> >::value_type neededMeshPair;
    for (const neededMeshPair& neededMeshes : _meshConfiguration->getNeededMeshes()){
      bool participantFound = false;
      for (const impl::PtrParticipant& participant : _participantConfiguration->getParticipants()){
        if(participant->getName()==neededMeshes.first){
          for (const std::string& neededMesh  : neededMeshes.second){
            bool meshFound = false;
            for (impl::MeshContext* meshContext : participant->usedMeshContexts()){
              if(meshContext->mesh->getName()==neededMesh){
                meshFound = true;
                break;
              }
            }
            PRECICE_CHECK(meshFound,
                  "The participant "<< neededMeshes.first <<
                  " needs to use the mesh " << neededMesh <<
                  " if he wants to use it in the coupling scheme.");
          }
          participantFound = true;
          break;
        }
      }
      PRECICE_ASSERT(participantFound);
    }

  }
}

int SolverInterfaceConfiguration:: getDimensions() const
{
  return _dimensions;
}

const PtrParticipantConfiguration &
SolverInterfaceConfiguration:: getParticipantConfiguration() const
{
  return _participantConfiguration;
}

}} // close namespaces


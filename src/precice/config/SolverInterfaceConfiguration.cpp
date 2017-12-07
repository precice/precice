#include "SolverInterfaceConfiguration.hpp"
#include "ParticipantConfiguration.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "xml/ValidatorEquals.hpp"
#include "xml/ValidatorOr.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"

namespace precice {
namespace config {

logging::Logger SolverInterfaceConfiguration:: _log("config::SolverInterfaceConfiguration");

SolverInterfaceConfiguration:: SolverInterfaceConfiguration
(
  utils::XMLTag& parent )
:
  TAG("solver-interface"),
  ATTR_DIMENSIONS("dimensions"),
  _dimensions(-1),
  _dataConfiguration(),
  _meshConfiguration(),
  _m2nConfiguration(),
  _participantConfiguration(),
  _couplingSchemeConfiguration()
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_ONCE);
  tag.setDocumentation("Configuration of simulation relevant features.");

  XMLAttribute<int> attrDimensions(ATTR_DIMENSIONS);
  doc = "Determines the spatial dimensionality of the configuration";
  attrDimensions.setDocumentation(doc);
  ValidatorEquals<int> validDim2(2);
  ValidatorEquals<int> validDim3(3);
  attrDimensions.setValidator(validDim2 || validDim3);
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
  utils::XMLTag& tag )
{
  TRACE();
  if (tag.getName() == TAG){
    _dimensions = tag.getIntAttributeValue(ATTR_DIMENSIONS);
    _dataConfiguration->setDimensions(_dimensions);
    _meshConfiguration->setDimensions(_dimensions);
    _participantConfiguration->setDimensions(_dimensions);
  }
  else {
    ERROR("Received callback from tag " << tag.getName());
  }
}

void SolverInterfaceConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  TRACE();
  if (tag.getName() == TAG){
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
            CHECK(meshFound,
                  "The participant "<< neededMeshes.first <<
                  " needs to use the mesh " << neededMesh <<
                  " if he wants to use it in the coupling scheme.");
          }
          participantFound = true;
          break;
        }
      }
      assertion(participantFound);
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


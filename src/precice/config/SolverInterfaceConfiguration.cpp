// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverInterfaceConfiguration.hpp"
#include "ParticipantConfiguration.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "spacetree/config/SpacetreeConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "geometry/Geometry.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/UncoupledScheme.hpp"
#include <limits>

namespace precice {
namespace config {

tarch::logging::Log SolverInterfaceConfiguration:: _log("precice::config::SolverInterfaceConfiguration");

SolverInterfaceConfiguration:: SolverInterfaceConfiguration
(
  utils::XMLTag& parent )
:
  TAG("solver-interface"),
  ATTR_DIMENSIONS("dimensions"),
  ATTR_GEOMETRY_MODE("geometry-mode"),
  ATTR_RESTART_MODE("restart-mode"),
  ATTR_SPACETREE_NAME("name"),
  _dimensions(-1),
  _geometryMode(false),
  _restartMode(false),
  //_participants(),
  //_indexAccessor(-1),
  _dataConfiguration(),
  _meshConfiguration(),
  _m2nConfiguration(),
  _geometryConfiguration(),
  _spacetreeConfiguration(),
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

  XMLAttribute<bool> attrGeometryMode(ATTR_GEOMETRY_MODE);
  doc = "By default it is assumed that preCICE is used to perform partitioned ";
  doc += "coupled simulations. If geometry-mode is activated, preCICE can be ";
  doc += "used by a single solver only.";
  attrGeometryMode.setDocumentation(doc);
  attrGeometryMode.setDefaultValue(false);
  tag.addAttribute(attrGeometryMode);

  XMLAttribute<bool> attrRestartMode(ATTR_RESTART_MODE);
  doc = "If restart-mode is activated, a formerly created simulation checkpoint ";
  doc += "is read at start of the simulation. The participating solvers have to ";
  doc += "write/read there own checkpoints of simulation data.";
  attrRestartMode.setDocumentation(doc);
  attrRestartMode.setDefaultValue(false);
  tag.addAttribute(attrRestartMode);

  _dataConfiguration = mesh::PtrDataConfiguration (
      new mesh::DataConfiguration(tag) );
  _meshConfiguration = mesh::PtrMeshConfiguration (
      new mesh::MeshConfiguration(tag, _dataConfiguration) );
  _m2nConfiguration = m2n::M2NConfiguration::SharedPointer (
      new m2n::M2NConfiguration(tag) );
  _geometryConfiguration = geometry::PtrGeometryConfiguration (
      new geometry::GeometryConfiguration(tag, _meshConfiguration) );
  _spacetreeConfiguration = spacetree::PtrSpacetreeConfiguration (
      new spacetree::SpacetreeConfiguration(tag) );
  _participantConfiguration = config::PtrParticipantConfiguration (
      new ParticipantConfiguration(tag, _meshConfiguration,
     _geometryConfiguration, _spacetreeConfiguration) );
  _couplingSchemeConfiguration = cplscheme::PtrCouplingSchemeConfiguration (
    new cplscheme::CouplingSchemeConfiguration(tag, _meshConfiguration,
    _m2nConfiguration) );

  parent.addSubtag(tag);
}

void SolverInterfaceConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getName());
  if (tag.getName() == TAG){
    _dimensions = tag.getIntAttributeValue(ATTR_DIMENSIONS);
    _geometryMode = tag.getBooleanAttributeValue(ATTR_GEOMETRY_MODE);
    _restartMode = tag.getBooleanAttributeValue(ATTR_RESTART_MODE);
    _dataConfiguration->setDimensions(_dimensions);
    _meshConfiguration->setDimensions(_dimensions);
    _geometryConfiguration->setDimensions(_dimensions);
    _spacetreeConfiguration->setDimensions(_dimensions);
    _participantConfiguration->setDimensions(_dimensions);
  }
  else {
    preciceError("xmlTagCallback()", "Received callback from tag " << tag.getName());
  }
}

void SolverInterfaceConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getName());
  if (tag.getName() == TAG){
    _meshConfiguration->setMeshSubIDs();
    if (_geometryMode ){
      assertion ( _participantConfiguration->getParticipants().size() == 1 );
      preciceCheck(not _participantConfiguration->getParticipants()[0]->useMaster(),
          "xmlEndTagCallback()", "In geometry mode, the usage of a master is not yet supported");
      std::string name = _participantConfiguration->getParticipants()[0]->getName();
      if ( not _couplingSchemeConfiguration->hasCouplingScheme(name)){
        double maxTime = cplscheme::CouplingScheme::UNDEFINED_TIME;
        int maxTimesteps = cplscheme::CouplingScheme::UNDEFINED_TIMESTEPS;
        int validDigits = 10;
        cplscheme::PtrCouplingScheme cplScheme (
            new cplscheme::UncoupledScheme(maxTime, maxTimesteps,
            validDigits, name) );
        _couplingSchemeConfiguration->addCouplingScheme ( cplScheme, name );
      }
    }
    else{

      //test if both participants do have the exchange meshes
      typedef std::map<std::string, std::vector<std::string> >::value_type neededMeshPair;
      foreach (const neededMeshPair& neededMeshes, _meshConfiguration->getNeededMeshes()){
        bool participantFound = false;
        foreach(const impl::PtrParticipant& participant, _participantConfiguration->getParticipants()){
          if(participant->getName()==neededMeshes.first){
            foreach(const std::string& neededMesh ,neededMeshes.second){
              bool meshFound = false;
              for (impl::MeshContext* meshContext : participant->usedMeshContexts()){
                if(meshContext->mesh->getName()==neededMesh){
                  meshFound = true;
                  break;
                }
              }
              preciceCheck(meshFound,"xmlEndTagCallback()",
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
}

int SolverInterfaceConfiguration:: getDimensions() const
{
  return _dimensions;
}

const spacetree::PtrSpacetreeConfiguration &
SolverInterfaceConfiguration:: getSpacetreeConfiguration() const
{
  return _spacetreeConfiguration;
}

const PtrParticipantConfiguration &
SolverInterfaceConfiguration:: getParticipantConfiguration() const
{
  return _participantConfiguration;
}

}} // close namespaces

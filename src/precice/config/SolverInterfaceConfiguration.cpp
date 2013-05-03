// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverInterfaceConfiguration.hpp"
#include "ParticipantConfiguration.hpp"
#include "precice/impl/Participant.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "spacetree/config/SpacetreeConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "geometry/Geometry.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/UncoupledCouplingScheme.hpp"
#include <limits>

namespace precice {
namespace config {

tarch::logging::Log SolverInterfaceConfiguration:: _log("precice::config::SolverInterfaceConfiguration");


//const std::string& SolverInterfaceConfiguration:: getTag()
//{
//  static std::string tag("solver-interface");
//  return tag;
//}

SolverInterfaceConfiguration:: SolverInterfaceConfiguration
(
  utils::XMLTag& parent )
:
  TAG("solver-interface"),
  ATTR_DIMENSIONS("dimensions"),
  ATTR_GEOMETRY_MODE("geometry-mode"),
  ATTR_RESTART_MODE("restart-mode"),
  ATTR_SPACETREE_NAME("name"),
  //_isValid(false),
  _geometryMode(false),
  _restartMode(false),
  _participants(),
  _dataConfiguration(),
  _meshConfiguration(),
  _comConfiguration(),
  _geometryConfiguration(),
  _spacetreeConfiguration(),
  _participantConfiguration(),
  _couplingSchemeConfiguration(),
  _indexAccessor(-1)
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
  _comConfiguration = com::PtrCommunicationConfiguration (
      new com::CommunicationConfiguration(tag) );
  _geometryConfiguration = geometry::PtrGeometryConfiguration (
      new geometry::GeometryConfiguration(tag, _meshConfiguration) );
  _spacetreeConfiguration = spacetree::PtrSpacetreeConfiguration (
      new spacetree::SpacetreeConfiguration(tag) );
  _participantConfiguration = config::PtrParticipantConfiguration (
      new ParticipantConfiguration(tag, _meshConfiguration,
     _geometryConfiguration, _spacetreeConfiguration) );
  _couplingSchemeConfiguration = cplscheme::PtrCouplingSchemeConfiguration (
    new cplscheme::CouplingSchemeConfiguration(tag, _meshConfiguration,
    _comConfiguration) );

  parent.addSubtag(tag);
}

//bool SolverInterfaceConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
  //_dimensions = readDimensions(xmlReader);
  //_geometryMode = readGeometryMode(xmlReader);
  //_restartMode = readRestartMode(xmlReader);
//
  //_dataConfiguration = mesh::PtrDataConfiguration (
  //    new mesh::DataConfiguration(_dimensions) );
  //_meshConfiguration = mesh::PtrMeshConfiguration (
  //    new mesh::MeshConfiguration(_dataConfiguration) );
  //_comConfiguration = com::PtrCommunicationConfiguration (
  //    new com::CommunicationConfiguration() );
  //_geometryConfiguration = geometry::PtrGeometryConfiguration (
  //    new geometry::GeometryConfiguration(_meshConfiguration, _dimensions) );
  //_spacetreeConfiguration = spacetree::PtrSpacetreeConfiguration (
  //    new spacetree::SpacetreeConfiguration(_dimensions) );
  //_participantConfiguration = config::PtrParticipantConfiguration (
  //    new ParticipantConfiguration(_dimensions, _meshConfiguration,
  //   _geometryConfiguration, _spacetreeConfiguration) );
  //_couplingSchemeConfiguration = cplscheme::PtrCouplingSchemeConfiguration (
  //  new cplscheme::CouplingSchemeConfiguration(_meshConfiguration,
  //  _comConfiguration) );
//
  //using namespace utils;
  //XMLTag tagSolverInterfaceImpl ( TAG, XMLTag::OCCUR_ONCE );
//
  //XMLTag subtagData ( mesh::DataConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
  //_tag.addSubtag ( subtagData );
//
  //XMLTag subtagCom (
  //    com::CommunicationConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
  //_tag.addSubtag ( subtagCom );
//
  //XMLTag::Occurrence occurrence = XMLTag::OCCUR_ARBITRARY;
  //if ( _geometryMode ){
  //  occurrence = XMLTag::OCCUR_NOT_OR_ONCE;
  //}
  //XMLTag subtagParticipant ( ParticipantConfiguration::TAG, occurrence );
  //_tag.addSubtag ( subtagParticipant );
//
  //XMLTag subtagMesh ( mesh::MeshConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
  //_tag.addSubtag ( subtagMesh );
//
  //XMLTag subtagGeometry (
  //    geometry::GeometryConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
  //_tag.addSubtag ( subtagGeometry );
//
  //XMLTag subtagSpacetree (
  //    spacetree::SpacetreeConfiguration::TAG, XMLTag::OCCUR_ARBITRARY );
  //XMLAttribute<std::string> name ( ATTR_SPACETREE_NAME );
  //subtagSpacetree.addAttribute ( name );
  //_tag.addSubtag ( subtagSpacetree );
//
  //if ( not _geometryMode ){
  //  XMLTag subtagCouplingScheme (
  //      cplscheme::CouplingSchemeConfiguration::TAG, XMLTag::OCCUR_ONCE );
  //  _tag.addSubtag ( subtagCouplingScheme );
  //}
//  _isValid = _tag.parse(xmlReader);
//
//  if (_isValid){
//    _meshConfiguration->setMeshSubIDs();
//  }
//  return _isValid;
//}

//bool SolverInterfaceConfiguration:: isValid() const
//{
//  return _isValid;
//}

void SolverInterfaceConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1 ( "xmlTagCallback()", tag.getName() );
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
//  if ( tag.getName() == mesh::DataConfiguration::TAG ) {
//    return _dataConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == mesh::MeshConfiguration::TAG ) {
//    return _meshConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == com::CommunicationConfiguration::TAG ) {
//    return _comConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == geometry::GeometryConfiguration::TAG ) {
//    return _geometryConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == spacetree::SpacetreeConfiguration::TAG ) {
//    return  _spacetreeConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == ParticipantConfiguration::TAG ) {
//    return _participantConfiguration->parseSubtag ( xmlReader );
//  }
//  else if ( tag.getName() == cplscheme::CouplingSchemeConfiguration::TAG ) {
//    return _couplingSchemeConfiguration->parseSubtag ( xmlReader );
//  }
}

void SolverInterfaceConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  if (tag.getName() == TAG){
    _meshConfiguration->setMeshSubIDs();
    if (_geometryMode ){
      assertion ( _participantConfiguration->getParticipants().size() == 1 );
      std::string name = _participantConfiguration->getParticipants()[0]->getName();
      if ( not _couplingSchemeConfiguration->hasCouplingScheme(name)){
        double maxTime = cplscheme::CouplingScheme::UNDEFINED_TIME;
        int maxTimesteps = cplscheme::CouplingScheme::UNDEFINED_TIMESTEPS;
        int validDigits = 10;
        cplscheme::PtrCouplingScheme cplScheme (
            new cplscheme::UncoupledCouplingScheme(maxTime, maxTimesteps,
            validDigits, name) );
        _couplingSchemeConfiguration->addCouplingScheme ( cplScheme, name );
      }
    }
  }
}

int SolverInterfaceConfiguration:: getDimensions() const
{
  //assertion ( _isValid );
  return _dimensions;
}

const spacetree::PtrSpacetreeConfiguration &
SolverInterfaceConfiguration:: getSpacetreeConfiguration () const
{
  return _spacetreeConfiguration;
}

const PtrParticipantConfiguration &
SolverInterfaceConfiguration:: getParticipantConfiguration () const
{
  return _participantConfiguration;
}

//int SolverInterfaceConfiguration:: readDimensions
//(
//  tarch::irr::io::IrrXMLReader* xmlReader )
//{
//  int attributeCount = xmlReader->getAttributeCount();
//  for ( int i=0; i < attributeCount; i++ ) {
//    if ( xmlReader->getAttributeName(i) == ATTR_DIMENSIONS ) {
//      int dimensions = xmlReader->getAttributeValueAsInt(i);
//      preciceCheck ( (dimensions == 2) || (dimensions == 3), "readDimensions()",
//          "Attribute \"dimensions\" of tag <" << TAG << "> has to be either"
//          << " 2 or 3!");
//      return dimensions;
//    }
//  }
//  preciceError ( "readDimensions()", "Attribute \"dimensions\" of tag <"
//      << TAG << "> is missing!" );
//}

//bool SolverInterfaceConfiguration:: readGeometryMode
//(
//   tarch::irr::io::IrrXMLReader * xmlReader )
//{
//  int attributeCount = xmlReader->getAttributeCount();
//  for ( int i=0; i < attributeCount; i++ ) {
//    if ( xmlReader->getAttributeName(i) == ATTR_GEOMETRY_MODE ) {
//      return xmlReader->getAttributeValueAsBool ( i );
//    }
//  }
//  return false;
//}

//bool SolverInterfaceConfiguration:: readRestartMode
//(
//  tarch::irr::io::IrrXMLReader * xmlReader )
//{
//  int attributeCount = xmlReader->getAttributeCount();
//  for ( int i=0; i < attributeCount; i++ ) {
//    if ( xmlReader->getAttributeName(i) == ATTR_RESTART_MODE ) {
//      return xmlReader->getAttributeValueAsBool ( i );
//    }
//  }
//  return false;
//}

}} // close namespaces

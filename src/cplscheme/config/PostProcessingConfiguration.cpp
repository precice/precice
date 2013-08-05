// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PostProcessingConfiguration.hpp"
#include "cplscheme/impl/ConstantRelaxationPostProcessing.hpp"
#include "cplscheme/impl/AitkenPostProcessing.hpp"
#include "cplscheme/impl/HierarchicalAitkenPostProcessing.hpp"
#include "cplscheme/impl/IQNILSPostProcessing.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log PostProcessingConfiguration::
      _log ( "precice::cplscheme::PostProcessingConfiguration" );

//const std::string & PostProcessingConfiguration:: getTag ()
//{
//  static std::string tag ( "post-processing" );
//  return tag;
//}

PostProcessingConfiguration:: PostProcessingConfiguration
(
  const mesh::PtrMeshConfiguration& meshConfig )
:
  TAG("post-processing"),
  TAG_RELAX("relaxation" ),
  TAG_INIT_RELAX("initial-relaxation"),
  TAG_MAX_USED_ITERATIONS("max-used-iterations"),
  TAG_TIMESTEPS_REUSED("timesteps-reused"),
  TAG_SINGULARITY_LIMIT("singularity-limit"),
  ATTR_DATA("data"),
  ATTR_MESH("mesh"),
  ATTR_VALUE("value"),
  VALUE_CONSTANT("constant"),
  VALUE_AITKEN ("aitken"),
  VALUE_HIERARCHICAL_AITKEN("hierarchical-aitken"),
  VALUE_IQNILS ("IQN-ILS"),
  //_isValid(false),
  _meshConfig(meshConfig),
  _postProcessing(),
  _config()
{
  assertion(meshConfig.get() != NULL);
}

void PostProcessingConfiguration:: connectTags(
    utils::XMLTag&                    parent){

  using namespace utils;

    XMLTag::Occurrence occ = XMLTag::OCCUR_NOT_OR_ONCE;
    std::list<XMLTag> tags;
    {
      XMLTag tag(*this, VALUE_CONSTANT, occ, TAG);
      addTypeSpecificSubtags(tag);
      tags.push_back(tag);
    }
    {
      XMLTag tag(*this, VALUE_AITKEN, occ, TAG);
      addTypeSpecificSubtags(tag);
      tags.push_back(tag);
    }
    {
      XMLTag tag(*this, VALUE_HIERARCHICAL_AITKEN, occ, TAG);
      addTypeSpecificSubtags(tag);
      tags.push_back(tag);
    }
    {
      XMLTag tag(*this, VALUE_IQNILS, occ, TAG);
      addTypeSpecificSubtags(tag);
      tags.push_back(tag);
    }

    XMLAttribute<std::string> attrData(ATTR_DATA);
    XMLAttribute<std::string> attrMesh(ATTR_MESH);

    foreach (XMLTag& tag, tags){
      tag.addAttribute(attrData);
      tag.addAttribute(attrMesh);
      parent.addSubtag(tag);
    }
}

//bool PostProcessingConfiguration:: parseSubtag
//(
//  tarch::irr::io::IrrXMLReader * xmlReader )
//{
//  preciceTrace("parseSubtag()" );
//  using namespace utils;
//
//  XMLTag tag(TAG, XMLTag::OCCUR_ONCE );
//
//  XMLAttribute<std::string> attrType(ATTR_TYPE );
//  ValidatorEquals<std::string> validConstant(VALUE_CONSTANT );
//  ValidatorEquals<std::string> validAitken(VALUE_AITKEN );
//  ValidatorEquals<std::string> validHierarchAitken(VALUE_HIERARCHICAL_AITKEN );
//  ValidatorEquals<std::string> validIQNILS(VALUE_IQNILS );
//  attrType.setValidator (
//    validConstant || validAitken || validHierarchAitken|| validIQNILS );
//  tag.addAttribute(attrType );
//
//  XMLAttribute<std::string> attrData(ATTR_DATA );
//  tag.addAttribute(attrData );
//
//  XMLAttribute<std::string> attrMesh(ATTR_MESH );
//  tag.addAttribute(attrMesh );
//
//  _isValid = tag.parse(xmlReader, *this );
//  _config = ConfigurationData ();
//  return _isValid;
//}

//bool PostProcessingConfiguration:: isValid() const
//{
//  return _isValid;
//}

impl::PtrPostProcessing PostProcessingConfiguration:: getPostProcessing()
{
  return _postProcessing;
}

void PostProcessingConfiguration:: xmlTagCallback
(
  utils::XMLTag& callingTag )
{
  preciceTrace1("xmlTagCallback()", callingTag.getFullName());
  if (callingTag.getNamespace() == TAG){
    std::string dataName = callingTag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = callingTag.getStringAttributeValue(ATTR_MESH);
    foreach(mesh::PtrMesh mesh, _meshConfig->meshes() ) {
      if(mesh->getName() == meshName ) {
        foreach(mesh::PtrData data, mesh->data() ) {
          if (dataName == data->getName()){
            assertion(_config.dataID == -1 );
            _config.dataID = data->getID();
          }
        }
      }
    }
    if (_config.dataID == -1){
      std::ostringstream stream;
      stream << "Data with name \"" << dataName << "\" associated to mesh \""
             << meshName << "\" not found on configuration of post-processing";
      throw stream.str();
    }
    //_config.type = callingTag.getName();
    //addTypeSpecificSubtags(callingTag );
  }
  else if (callingTag.getName() == TAG_RELAX){
    _config.relaxationFactor = callingTag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (callingTag.getName() == TAG_INIT_RELAX){
    _config.relaxationFactor = callingTag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (callingTag.getName() == TAG_MAX_USED_ITERATIONS){
    _config.maxIterationsUsed = callingTag.getIntAttributeValue(ATTR_VALUE);
  }
  else if (callingTag.getName() == TAG_TIMESTEPS_REUSED){
    _config.timestepsReused = callingTag.getIntAttributeValue(ATTR_VALUE);
  }
  else if (callingTag.getName() == TAG_SINGULARITY_LIMIT){
    _config.singularityLimit = callingTag.getDoubleAttributeValue(ATTR_VALUE);
  }
}

void PostProcessingConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& callingTag )
{
  preciceTrace1("xmlEndTagCallback()", callingTag.getName());
  if (callingTag.getNamespace() == TAG){
    if (callingTag.getName() == VALUE_CONSTANT){
      _postProcessing = impl::PtrPostProcessing (
          new impl::ConstantRelaxationPostProcessing (
          _config.relaxationFactor, _config.dataID) );
    }
    else if (callingTag.getName() == VALUE_AITKEN){
      _postProcessing = impl::PtrPostProcessing (
          new impl::AitkenPostProcessing(
          _config.relaxationFactor, _config.dataID) );
    }
    else if (callingTag.getName() == VALUE_HIERARCHICAL_AITKEN){
      _postProcessing = impl::PtrPostProcessing (
          new impl::HierarchicalAitkenPostProcessing (
          _config.relaxationFactor, _config.dataID) );
    }
    else if (callingTag.getName() == VALUE_IQNILS){
      _postProcessing = impl::PtrPostProcessing (
          new impl::IQNILSPostProcessing(
          _config.relaxationFactor, _config.maxIterationsUsed,
          _config.timestepsReused, _config.singularityLimit,
          _config.dataID) );
    }
    else {
      assertion(false );
    }
  }
}

void PostProcessingConfiguration:: clear()
{
  _config = ConfigurationData();
  _postProcessing = impl::PtrPostProcessing();
}

void PostProcessingConfiguration:: addTypeSpecificSubtags
(
  utils::XMLTag& tag )
{
  using namespace utils;
  if(tag.getName() == VALUE_CONSTANT ) {
    XMLTag tagRelax(*this, TAG_RELAX, XMLTag::OCCUR_ONCE );
    XMLAttribute<double> attrValue(ATTR_VALUE );
    tagRelax.addAttribute(attrValue );
    tag.addSubtag(tagRelax );
  }
  else if (tag.getName() == VALUE_AITKEN){
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE );
    XMLAttribute<double> attrValue(ATTR_VALUE );
    tagInitRelax.addAttribute(attrValue );
    tag.addSubtag(tagInitRelax );
  }
  else if (tag.getName() == VALUE_HIERARCHICAL_AITKEN){
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE );
    XMLAttribute<double> attrValue(ATTR_VALUE );
    tagInitRelax.addAttribute(attrValue );
    tag.addSubtag(tagInitRelax );
  }
  else if (tag.getName() == VALUE_IQNILS){
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE );
    XMLAttribute<double> attrDoubleValue(ATTR_VALUE);
    tagInitRelax.addAttribute(attrDoubleValue);
    tag.addSubtag(tagInitRelax);

    XMLTag tagMaxUsedIter(*this, TAG_MAX_USED_ITERATIONS, XMLTag::OCCUR_ONCE );
    XMLAttribute<int> attrIntValue(ATTR_VALUE );
    tagMaxUsedIter.addAttribute(attrIntValue );
    tag.addSubtag(tagMaxUsedIter );

    XMLTag tagTimestepsReused(*this, TAG_TIMESTEPS_REUSED, XMLTag::OCCUR_ONCE );
    tagTimestepsReused.addAttribute(attrIntValue );
    tag.addSubtag(tagTimestepsReused );

    XMLTag tagSingularityLimit(*this, TAG_SINGULARITY_LIMIT, XMLTag::OCCUR_ONCE );
    tagSingularityLimit.addAttribute(attrDoubleValue );
    tag.addSubtag(tagSingularityLimit );
  }
  else {
    preciceError("addTypeSpecificSubtag()", "Post-processing of type \""
                 << tag.getName() << "\" is unknown!" );
  }
}

}} // namespace precice, cplscheme

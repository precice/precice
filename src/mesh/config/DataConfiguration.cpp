// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "DataConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace mesh {

tarch::logging::Log DataConfiguration:: _log("precice::mesh::DataConfiguration");

//const std::string& DataConfiguration:: getTag()
//{
//  static std::string tag("data");
//  return tag;
//}

DataConfiguration:: DataConfiguration
(
  utils::XMLTag& parent )
:
  TAG("data"),
  ATTR_NAME("name"),
  //ATTR_TYPE("type"),
  VALUE_VECTOR("vector"),
  VALUE_SCALAR("scalar"),
  _dimensions(0),
  //_isValid(false),
  _data(),
  _indexLastConfigured(-1)
{
  using namespace utils;
  std::string doc;
  XMLTag tagScalar(*this, VALUE_SCALAR, XMLTag::OCCUR_ARBITRARY, TAG);
  doc = "Defines a scalar data set to be assigned to meshes.";
  tagScalar.setDocumentation(doc);
  XMLTag tagVector(*this, VALUE_VECTOR, XMLTag::OCCUR_ARBITRARY, TAG);
  doc = "Defines a vector data set to be assigned to meshes. The number of ";
  doc += "components of each data entry depends on the spatial dimensions set ";
  doc += "in tag <solver-interface>.";
  tagVector.setDocumentation(doc);

  XMLAttribute<std::string> attrName(ATTR_NAME);
  doc = "Unique name for the data set.";
  attrName.setDocumentation(doc);
  tagScalar.addAttribute(attrName);
  tagVector.addAttribute(attrName);

//  XMLAttribute<std::string> attrType(ATTR_TYPE);
//  ValidatorEquals<std::string> validVectorType(VALUE_VECTOR);
//  ValidatorEquals<std::string> validFloatingType(VALUE_SCALAR);
//  attrType.setValidator(validVectorType || validFloatingType);
//  tag.addAttribute(attrType);

  parent.addSubtag(tagScalar);
  parent.addSubtag(tagVector);
}

void DataConfiguration:: setDimensions
(
  int dimensions )
{
  preciceTrace1("setDimensions()", dimensions);
  assertion1((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

//bool DataConfiguration:: parseSubtag
//(
//  tarch::irr::io::IrrXMLReader* xmlReader )
//{
  //utils::XMLTag tagData ( TAG, utils::XMLTag::OCCUR_ONCE );
//
//  utils::XMLAttribute<std::string> attrName ( ATTR_NAME );
//  tagData.addAttribute ( attrName );
//
//  utils::XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  utils::ValidatorEquals<std::string> validVectorType ( VALUE_VECTOR );
//  utils::ValidatorEquals<std::string> validFloatingType ( VALUE_SCALAR );
//  attrType.setValidator ( validVectorType || validFloatingType );
//  tagData.addAttribute ( attrType );
//
//  _isValid = _tag.parse(xmlReader);
//  return _isValid;
//}

//bool DataConfiguration:: isValid() const
//{
//   return _isValid;
//}

const std::vector<DataConfiguration::ConfiguredData>&
DataConfiguration:: data() const
{
   return _data;
}

DataConfiguration::ConfiguredData DataConfiguration:: getRecentlyConfiguredData() const
{
  assertion(_data.size() > 0);
  assertion(_indexLastConfigured >= 0);
  assertion(_indexLastConfigured < (int)_data.size());
  return _data[_indexLastConfigured];
}

void DataConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  if (tag.getNamespace() == TAG){
    assertion(_dimensions != 0);
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    std::string typeName = tag.getName();
    int dataDimensions = getDataDimensions(typeName);
    addData(name, dataDimensions);
  }
  else {
    preciceError("xmlTagCallback()", "Received callback from tag " << tag.getName());
  }
}

void DataConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
}

void DataConfiguration:: addData
(
  const std::string& name,
  int                dataDimensions )
{
  ConfiguredData data(name, dataDimensions);

  // Check, if data with same name has been added already
  for ( size_t i=0; i < _data.size(); i++ ) {
    preciceCheck ( _data[i].name != data.name, "addData()",
                   "Data \"" << data.name << "\" uses non-unique name!" );
  }
  _data.push_back ( data );
}

int DataConfiguration:: getDataDimensions
(
  const std::string& typeName ) const
{
  if (typeName == VALUE_VECTOR) {
    return _dimensions;
  }
  else if (typeName == VALUE_SCALAR) {
    return 1;
  }
  preciceError ( "getDataDimension()", "Unknown data type!" );
}

//int DataConfiguration:: getDimensions() const
//{
//  return _dimensions;
//}


}} // namespace precice, mesh

#include "DataConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/PropertyContainer.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/ValidatorEquals.hpp"
#include "xml/ValidatorOr.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace mesh {

logging::Logger DataConfiguration:: _log("mesh::DataConfiguration");


DataConfiguration:: DataConfiguration
(
  xml::XMLTag& parent )
:
  TAG("data"),
  ATTR_NAME("name"),
  VALUE_VECTOR("vector"),
  VALUE_SCALAR("scalar"),
  _dimensions(0),
  _data(),
  _indexLastConfigured(-1)
{
  using namespace xml;
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

  parent.addSubtag(tagScalar);
  parent.addSubtag(tagVector);
}

void DataConfiguration:: setDimensions
(
  int dimensions )
{
  TRACE(dimensions);
  assertion((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

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
  xml::XMLTag& tag )
{
  if (tag.getNamespace() == TAG){
    assertion(_dimensions != 0);
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    std::string typeName = tag.getName();
    int dataDimensions = getDataDimensions(typeName);
    addData(name, dataDimensions);
  }
  else {
    ERROR("Received callback from tag " << tag.getName());
  }
}

void DataConfiguration:: xmlEndTagCallback
(
  xml::XMLTag& tag )
{
}

void DataConfiguration:: addData
(
  const std::string& name,
  int                dataDimensions )
{
  ConfiguredData data(name, dataDimensions);

  // Check, if data with same name has been added already
  for (auto & elem : _data) {
    CHECK ( elem.name != data.name, "Data \"" << data.name << "\" uses non-unique name!" );
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
  ERROR("Unknown data type!" );
}



}} // namespace precice, mesh

#include "DataConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/PropertyContainer.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace mesh {

DataConfiguration:: DataConfiguration(xml::XMLTag& parent)
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

 auto attrName = XMLAttribute<std::string>(ATTR_NAME)
      .setDocumentation("Unique name for the data set.");
  tagScalar.addAttribute(attrName);
  tagVector.addAttribute(attrName);

  parent.addSubtag(tagScalar);
  parent.addSubtag(tagVector);
}

void DataConfiguration:: setDimensions
(
  int dimensions )
{
  P_TRACE(dimensions);
  P_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

const std::vector<DataConfiguration::ConfiguredData>&
DataConfiguration:: data() const
{
   return _data;
}

DataConfiguration::ConfiguredData DataConfiguration:: getRecentlyConfiguredData() const
{
  P_ASSERT(_data.size() > 0);
  P_ASSERT(_indexLastConfigured >= 0);
  P_ASSERT(_indexLastConfigured < (int)_data.size());
  return _data[_indexLastConfigured];
}

void DataConfiguration:: xmlTagCallback
(
  xml::XMLTag& tag )
{
  if (tag.getNamespace() == TAG){
    P_ASSERT(_dimensions != 0);
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    std::string typeName = tag.getName();
    int dataDimensions = getDataDimensions(typeName);
    addData(name, dataDimensions);
  }
  else {
    P_ERROR("Received callback from tag " << tag.getName());
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
    P_CHECK( elem.name != data.name, "Data \"" << data.name << "\" uses non-unique name!" );
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
  P_ERROR("Unknown data type!" );
}



}} // namespace precice, mesh

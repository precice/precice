#include "DataConfiguration.hpp"
#include <ostream>
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::mesh {

DataConfiguration::DataConfiguration(xml::XMLTag &parent)
{
  using namespace xml;

  auto attrName = XMLAttribute<std::string>(ATTR_NAME)
                      .setDocumentation("Unique name for the data set.");

  auto attrDegree = makeXMLAttribute(ATTR_DEGREE, time::Time::DEFAULT_WAVEFORM_DEGREE);
  attrDegree.setDocumentation("Polynomial degree of waveform that is used for time interpolation.");

  XMLTag tagScalar(*this, VALUE_SCALAR, XMLTag::OCCUR_ARBITRARY, TAG);
  tagScalar.setDocumentation("Defines a scalar data set to be assigned to meshes.");
  tagScalar.addAttribute(attrName);
  tagScalar.addAttribute(attrDegree);
  parent.addSubtag(tagScalar);

  XMLTag tagVector(*this, VALUE_VECTOR, XMLTag::OCCUR_ARBITRARY, TAG);
  tagVector.setDocumentation("Defines a vector data set to be assigned to meshes. The number of "
                             "components of each data entry depends on the spatial dimensions of the mesh.");
  tagVector.addAttribute(attrName);
  tagVector.addAttribute(attrDegree);
  parent.addSubtag(tagVector);
}

const std::vector<DataConfiguration::ConfiguredData> &
DataConfiguration::data() const
{
  return _data;
}

DataConfiguration::ConfiguredData DataConfiguration::getRecentlyConfiguredData() const
{
  PRECICE_ASSERT(_data.size() > 0);
  PRECICE_ASSERT(_indexLastConfigured >= 0);
  PRECICE_ASSERT(_indexLastConfigured < (int) _data.size());
  return _data[_indexLastConfigured];
}

void DataConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  if (tag.getNamespace() == TAG) {
    const std::string &name = tag.getStringAttributeValue(ATTR_NAME);

    Data::typeName typeName;
    if (tag.getName() == "scalar") {
      typeName = Data::typeName::SCALAR;
    } else if (tag.getName() == "vector") {
      typeName = Data::typeName::VECTOR;
    } else {
      PRECICE_ERROR("You configured data with name=\"{}\" to be of type \"{}\", but this type is unknown. Known types are \"scalar\" and \"vector\".", name, tag.getName());
    };

    const int waveformDegree = tag.getIntAttributeValue(ATTR_DEGREE);
    PRECICE_CHECK(!(waveformDegree < time::Time::MIN_WAVEFORM_DEGREE),
                  "You tried to configure the data with name \"{}\" to use the waveform-degree=\"{}\", but the degree must be at least \"{}\".", name, waveformDegree, time::Time::MIN_WAVEFORM_DEGREE);
    addData(name, typeName, waveformDegree);
  } else {
    PRECICE_ASSERT(false, "Received callback from an unknown tag.", tag.getName());
  }
}

void DataConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
}

void DataConfiguration::addData(
    const std::string &  name,
    const Data::typeName typeName,
    int                  waveformDegree)
{
  // Check if data with same name has been added already
  for (auto &elem : _data) {
    PRECICE_CHECK(elem.name != name,
                  "Data \"{0}\" has already been defined. Please rename or remove one of the data tags with name=\"{0}\".",
                  name);
  }

  _data.emplace_back(name, typeName, waveformDegree);
}

} // namespace precice::mesh

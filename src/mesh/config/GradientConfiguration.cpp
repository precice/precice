#include "GradientConfiguration.hpp"
#include <ostream>
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace mesh {

GradientConfiguration::GradientConfiguration(xml::XMLTag &parent)
{
  using namespace xml;
  std::string doc;
  XMLTag      tagScalar(*this, VALUE_SCALAR, XMLTag::OCCUR_ARBITRARY, TAG);
  doc = "Defines a scalar gradient data set to be assigned to meshes.";
  tagScalar.setDocumentation(doc);
  XMLTag tagVector(*this, VALUE_VECTOR, XMLTag::OCCUR_ARBITRARY, TAG);
  doc = "Defines a vector gradient data set to be assigned to meshes. The number of ";
  doc += "components of each gradient entry depends on the spatial dimensions set ";
  doc += "in tag <solver-interface>.";
  tagVector.setDocumentation(doc);

  auto attrName = XMLAttribute<std::string>(ATTR_NAME)
                      .setDocumentation("Unique name for the gradient data set.");
  tagScalar.addAttribute(attrName);
  tagVector.addAttribute(attrName);

  parent.addSubtag(tagScalar);
  parent.addSubtag(tagVector);
}

void GradientConfiguration::setDimensions(
    int dimensions)
{
  PRECICE_TRACE(dimensions);
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

const std::vector<GradientConfiguration::ConfiguredGradient> &
GradientConfiguration::gradients() const
{
  return _gradients;
}

GradientConfiguration::ConfiguredGradient GradientConfiguration::getRecentlyConfiguredGradient() const
{
  PRECICE_ASSERT(_gradients.size() > 0);
  PRECICE_ASSERT(_indexLastConfigured >= 0);
  PRECICE_ASSERT(_indexLastConfigured < (int) _gradients.size());
  return _gradients[_indexLastConfigured];
}

void GradientConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  if (tag.getNamespace() == TAG) {
    PRECICE_ASSERT(_dimensions != 0);
    std::string name           = tag.getStringAttributeValue(ATTR_NAME);
    std::string typeName       = tag.getName();
    int         dataDimensions = getDataDimensions(typeName);
    addGradient(name, dataDimensions);
  } else {
    PRECICE_ASSERT(false, "Received callback from unknown tag " << tag.getName());
  }
}

void GradientConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
}

void GradientConfiguration::addGradient(
    const std::string &name,
    int                dataDimensions)
{
  ConfiguredGradient gradient(name, dataDimensions);

  // Check if gradient with same name has been added already
  for (auto &elem : _gradients) {
    PRECICE_CHECK(elem.name != gradient.name, "Gradient \"" << gradient.name << "\" has already been defined. Please rename or remove one of the gradient tags with name=\"" << gradient.name << "\".");
  }
  _gradients.push_back(gradient);
}

int GradientConfiguration::getDataDimensions(
    const std::string &typeName) const
{
  if (typeName == VALUE_VECTOR) {
    return _dimensions;
  } else if (typeName == VALUE_SCALAR) {
    return 1;
  }
  PRECICE_ASSERT(false, "Unknown gradient data type \"" << typeName << "\". Known data types: " << VALUE_SCALAR << ", " << VALUE_VECTOR << ".");
}

} // namespace mesh
} // namespace precice

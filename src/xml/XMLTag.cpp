#include "xml/XMLTag.hpp"
#include <Eigen/Core>
#include <ostream>
#include <utility>
#include "logging/LogMacros.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"

namespace precice {
namespace xml {

XMLTag::XMLTag(
    Listener &  listener,
    std::string tagName,
    Occurrence  occurrence,
    std::string xmlNamespace)
    : _listener(listener),
      _name(std::move(tagName)),
      _namespace(std::move(xmlNamespace)),
      _occurrence(occurrence)
{
  if (not _namespace.empty()) {
    _fullName = _namespace + ":" + _name;
  } else {
    _fullName = _name;
  }
}

XMLTag &XMLTag::setDocumentation(const std::string &documentation)
{
  _doc = documentation;
  return *this;
}

XMLTag &XMLTag::addNamespace(const std::string &namespaceName)
{
  _namespaces.push_back(namespaceName);
  return *this;
}

XMLTag &XMLTag::addSubtag(const XMLTag &tag)
{
  PRECICE_TRACE(tag._fullName);
  PRECICE_ASSERT(tag._name != std::string(""));
  if (not tag._namespace.empty()) {
    _configuredNamespaces[tag._namespace] = false;
  }

  _subtags.push_back(std::make_shared<XMLTag>(tag));
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<double> &attribute)
{
  PRECICE_TRACE(attribute.getName());
  PRECICE_ASSERT(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _doubleAttributes.insert(std::pair<std::string, XMLAttribute<double>>(attribute.getName(), attribute));
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<int> &attribute)
{
  PRECICE_TRACE(attribute.getName());
  PRECICE_ASSERT(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _intAttributes.insert(std::pair<std::string, XMLAttribute<int>>(attribute.getName(), attribute));
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<std::string> &attribute)
{
  PRECICE_TRACE(attribute.getName());
  PRECICE_ASSERT(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _stringAttributes.insert(std::pair<std::string, XMLAttribute<std::string>>(attribute.getName(), attribute));
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<bool> &attribute)
{
  PRECICE_TRACE(attribute.getName());
  PRECICE_ASSERT(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _booleanAttributes.insert(std::pair<std::string, XMLAttribute<bool>>(attribute.getName(), attribute));
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<Eigen::VectorXd> &attribute)
{
  PRECICE_TRACE(attribute.getName());
  PRECICE_ASSERT(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _eigenVectorXdAttributes.insert(
      std::pair<std::string, XMLAttribute<Eigen::VectorXd>>(attribute.getName(), attribute));
  return *this;
}

bool XMLTag::hasAttribute(const std::string &attributeName)
{
  return utils::contained(attributeName, _attributes);
}

double XMLTag::getDoubleAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<double>>::const_iterator iter;
  iter = _doubleAttributes.find(name);
  PRECICE_ASSERT(iter != _doubleAttributes.end());
  return iter->second.getValue();
}

int XMLTag::getIntAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<int>>::const_iterator iter;
  iter = _intAttributes.find(name);
  PRECICE_ASSERT(iter != _intAttributes.end());
  return iter->second.getValue();
}

const std::string &XMLTag::getStringAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<std::string>>::const_iterator iter;
  iter = _stringAttributes.find(name);
  PRECICE_ASSERT(iter != _stringAttributes.end(), name);
  return iter->second.getValue();
}

bool XMLTag::getBooleanAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<bool>>::const_iterator iter;
  iter = _booleanAttributes.find(name);
  PRECICE_ASSERT(iter != _booleanAttributes.end());
  return iter->second.getValue();
}

Eigen::VectorXd XMLTag::getEigenVectorXdAttributeValue(const std::string &name, int dimensions) const
{
  PRECICE_TRACE(name, dimensions);
  // std::map<std::string, XMLAttribute<utils::DynVector> >::const_iterator iter;
  auto iter = _eigenVectorXdAttributes.find(name);
  PRECICE_ASSERT(iter != _eigenVectorXdAttributes.end());
  const auto size = iter->second.getValue().size();
  PRECICE_CHECK(size == dimensions,
                "Vector attribute \"" << name << "\" of tag <" << getFullName() << "> is "
                                      << size << "D, which does not match the dimension of the " << dimensions << "D solver-interface.");

  // Read only first "dimensions" components of the parsed vector values
  Eigen::VectorXd        result(dimensions);
  const Eigen::VectorXd &parsed = iter->second.getValue();
  for (int i = 0; i < dimensions; i++) {
    result[i] = parsed[i];
  }
  PRECICE_DEBUG("Returning value = " << result);
  return result;
}

void XMLTag::readAttributes(const std::map<std::string, std::string> &aAttributes)
{
  PRECICE_TRACE();

  for (auto &element : aAttributes) {
    auto name = element.first;

    if (not utils::contained(name, _attributes)) {
      PRECICE_ERROR("Tag <" << _name << "> contains an unknown attribute named \"" << name << "\".");
    }
  }

  for (auto &pair : _doubleAttributes) {
    pair.second.readValue(aAttributes);
  }

  for (auto &pair : _intAttributes) {
    pair.second.readValue(aAttributes);
  }

  for (auto &pair : _stringAttributes) {
    pair.second.readValue(aAttributes);
  }

  for (auto &pair : _booleanAttributes) {
    pair.second.readValue(aAttributes);
  }

  for (auto &pair : _eigenVectorXdAttributes) {
    pair.second.readValue(aAttributes);
  }
}

/*void XMLTag:: readAttributes
(
  XMLReader* xmlReader )
{
  PRECICE_TRACE();
//  using utils::contained;
//  std::set<std::string> readNames;
  for (int i=0; i < xmlReader->getAttributeCount(); i++){
    std::string name = xmlReader->getAttributeName(i);
    if (not utils::contained(name, _attributes)){
      std::string error = "Wrong attribute \"" + name + "\"";
      throw std::runtime_error{error};
    }
//    else if (contained(name, _doubleAttributes)){
//      XMLAttribute<double>& attr = _doubleAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _intAttributes)){
//      XMLAttribute<int>& attr = _intAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _stringAttributes)){
//      XMLAttribute<std::string>& attr = _stringAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _booleanAttributes)){
//      XMLAttribute<bool>& attr = _booleanAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _vector2DAttributes)){
//      XMLAttribute<Vector2D>& attr = _vector2DAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _vector3DAttributes)){
//      XMLAttribute<Vector3D>& attr = _vector3DAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else if (contained(name, _dynVectorAttributes)){
//      XMLAttribute<DynVector>& attr = _dynVectorAttributes[name];
//      attr.readValue(xmlReader);
//    }
//    else {
//      throw std::runtime_error{"Internal error in readAttributes"};
//    }
//    readNames.insert(name);
  }

//  // Check if all attributes are read
//  for (const std::string& name : _attributes){
//    if (not contained(name, readNames)){
//
//      std::ostringstream stream;
//      stream << "Attribute \"" << name << "\" is not defined";
//      throw std::runtime_error{stream.str()};
//    }
//  }

  for (auto & pair : _doubleAttributes){
     pair.second.readValue(xmlReader);
  }

  for ( auto & pair : _intAttributes){
    pair.second.readValue(xmlReader);
  }

  for ( auto & pair : _stringAttributes ){
    pair.second.readValue(xmlReader);
  }

  for ( auto & pair : _booleanAttributes ){
    pair.second.readValue(xmlReader);
  }

  for ( auto & pair : _eigenVectorXdAttributes ){
    pair.second.readValue(xmlReader);
  }
}*/

void XMLTag::areAllSubtagsConfigured() const
{
  for (auto tag : _subtags) {
    std::string ns         = tag->_namespace;
    bool        configured = tag->isConfigured();

    bool occurOnce       = tag->getOccurrence() == OCCUR_ONCE;
    bool occurOnceOrMore = tag->getOccurrence() == OCCUR_ONCE_OR_MORE;

    if (not ns.empty()) {
      PRECICE_ASSERT(utils::contained(ns, _configuredNamespaces));
      configured |= _configuredNamespaces.find(ns)->second;
    }

    if ((not configured) && (occurOnce || occurOnceOrMore)) {

      if (tag->getNamespace().empty()) {
        PRECICE_ERROR("Tag <" << tag->getName() << "> was not found but is required to occur at least once.");
      } else {
        PRECICE_ERROR("Tag <" << tag->getNamespace() << ":... > was not found but is required to occur at least once.");
      }
    }
  }
}

void XMLTag::resetAttributes()
{
  _configured = false;

  for (auto &pair : _configuredNamespaces) {
    pair.second = false;
  }

  for (auto &pair : _doubleAttributes) {
    pair.second.setRead(false);
  }

  for (auto &pair : _intAttributes) {
    pair.second.setRead(false);
  }

  for (auto &pair : _stringAttributes) {
    pair.second.setRead(false);
  }

  for (auto &pair : _booleanAttributes) {
    pair.second.setRead(false);
  }

  for (auto &pair : _eigenVectorXdAttributes) {
    pair.second.setRead(false);
  }

  for (auto &tag : _subtags) {
    tag->_configured = false;
    tag->resetAttributes();
  }
}

void XMLTag::clear()
{
  _doubleAttributes.clear();
  _intAttributes.clear();
  _stringAttributes.clear();
  _booleanAttributes.clear();
  _subtags.clear();
}

//NoPListener& getNoPListener()
//{
//  static NoPListener listener;
//  return listener;
//}

XMLTag getRootTag()
{
  static NoPListener listener;
  return XMLTag(listener, "configuration", XMLTag::OCCUR_ONCE);
}

void configure(
    XMLTag &                                  tag,
    const precice::xml::ConfigurationContext &context,
    const std::string &                       configurationFilename)
{
  logging::Logger _log("xml");
  PRECICE_TRACE(tag.getFullName(), configurationFilename);

  NoPListener nopListener;
  XMLTag      root(nopListener, "", XMLTag::OCCUR_ONCE);

  precice::xml::ConfigParser p(configurationFilename, context, std::make_shared<XMLTag>(tag));

  root.addSubtag(tag);
}

std::string XMLTag::getOccurrenceString(XMLTag::Occurrence occurrence)
{
  if (occurrence == XMLTag::OCCUR_ARBITRARY) {
    return std::string("0..*");
  } else if (occurrence == XMLTag::OCCUR_NOT_OR_ONCE) {
    return std::string("0..1");
  } else if (occurrence == XMLTag::OCCUR_ONCE) {
    return std::string("1");
  } else if (occurrence == XMLTag::OCCUR_ONCE_OR_MORE) {
    return std::string("1..*");
  }
  return "";
}
} // namespace xml
} // namespace precice

//std::ostream& operator<<
//(
//  std::ostream&                 os,
//  const precice::xml::XMLTag& tag )
//{
//  os << tag.printDocumentation(80, 0);
//  return os;
//}

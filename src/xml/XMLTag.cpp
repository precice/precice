#include "xml/XMLTag.hpp"
#include <Eigen/Core>
#include <utility>
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"

namespace precice::xml {

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

XMLTag &XMLTag::setDocumentation(std::string_view documentation)
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
  const auto &name = attribute.getName();
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributes.emplace_back(attribute);
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<int> &attribute)
{
  const auto &name = attribute.getName();
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributes.emplace_back(attribute);
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<std::string> &attribute)
{
  const auto &name = attribute.getName();
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributes.emplace_back(attribute);
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<bool> &attribute)
{
  const auto &name = attribute.getName();
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributes.emplace_back(attribute);
  return *this;
}

XMLTag &XMLTag::addAttribute(const XMLAttribute<Eigen::VectorXd> &attribute)
{
  const auto &name = attribute.getName();
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributes.emplace_back(attribute);
  return *this;
}

void XMLTag::addAttributeHint(std::string name, std::string message)
{
  PRECICE_TRACE(name);
  PRECICE_ASSERT(!hasAttribute(name) && _attributeHints.count(name) == 0);
  _attributeHints.emplace(std::move(name), std::move(message));
}

namespace {
auto findAttribute(const XMLTag::Attributes &attributes, const std::string &name)
{
  return std::find_if(attributes.begin(), attributes.end(), [&name](const auto &attribute) { return getName(attribute) == name; });
}
} // namespace

bool XMLTag::hasAttribute(const std::string &attributeName) const
{
  return findAttribute(_attributes, attributeName) != _attributes.end();
}

double XMLTag::getDoubleAttributeValue(const std::string &name, std::optional<double> default_value) const
{
  PRECICE_TRACE(name);
  if (auto iter = findAttribute(_attributes, name);
      iter != _attributes.end()) {
    PRECICE_ASSERT(std::holds_alternative<XMLAttribute<double>>(*iter));
    return std::get<XMLAttribute<double>>(*iter).getValue();
  }
  if (default_value) {
    return default_value.value();
  }
  PRECICE_UNREACHABLE("The XMLAttribute doesn't exist, check its default.");
}

int XMLTag::getIntAttributeValue(const std::string &name, std::optional<int> default_value) const
{
  PRECICE_TRACE(name);
  if (auto iter = findAttribute(_attributes, name);
      iter != _attributes.end()) {
    PRECICE_ASSERT(std::holds_alternative<XMLAttribute<int>>(*iter));
    return std::get<XMLAttribute<int>>(*iter).getValue();
  }
  if (default_value) {
    return default_value.value();
  }
  PRECICE_UNREACHABLE("The XMLAttribute doesn't exist, check its default.");
}

std::string XMLTag::getStringAttributeValue(const std::string &name, std::optional<std::string> default_value) const
{
  PRECICE_TRACE(name);
  if (auto iter = findAttribute(_attributes, name);
      iter != _attributes.end()) {
    PRECICE_ASSERT(std::holds_alternative<XMLAttribute<std::string>>(*iter));
    return std::get<XMLAttribute<std::string>>(*iter).getValue();
  }
  if (default_value) {
    return default_value.value();
  }
  PRECICE_UNREACHABLE("The XMLAttribute doesn't exist, check its default.");
}

bool XMLTag::getBooleanAttributeValue(const std::string &name, std::optional<bool> default_value) const
{
  PRECICE_TRACE(name);
  if (auto iter = findAttribute(_attributes, name);
      iter != _attributes.end()) {
    PRECICE_ASSERT(std::holds_alternative<XMLAttribute<bool>>(*iter));
    return std::get<XMLAttribute<bool>>(*iter).getValue();
  }
  if (default_value) {
    return default_value.value();
  }
  PRECICE_UNREACHABLE("The XMLAttribute doesn't exist, check its default.");
}

Eigen::VectorXd XMLTag::getEigenVectorXdAttributeValue(const std::string &name) const
{
  PRECICE_TRACE(name);
  if (auto iter = findAttribute(_attributes, name);
      iter != _attributes.end()) {
    PRECICE_ASSERT(std::holds_alternative<XMLAttribute<Eigen::VectorXd>>(*iter));
    return std::get<XMLAttribute<Eigen::VectorXd>>(*iter).getValue();
  }
  PRECICE_UNREACHABLE("The XMLAttribute doesn't exist, check its default.");
}

void XMLTag::readAttributes(const std::map<std::string, std::string> &aAttributes)
{
  PRECICE_TRACE();

  // Check for unexpected attributes and hints
  for (const auto &element : aAttributes) {
    const auto &name = element.first;
    if (hasAttribute(name)) {
      continue;
    }

    // check existing hints
    if (auto pos = _attributeHints.find(name);
        pos != _attributeHints.end()) {
      PRECICE_ERROR("The tag <{}> in the configuration contains the attribute \"{}\". {}", _fullName, name, pos->second);
    }

    auto expected = getAttributeNames();
    auto matches  = utils::computeMatches(name, expected);
    if (!matches.empty() && matches.front().distance < 3) {
      matches.erase(std::remove_if(matches.begin(), matches.end(), [](auto &m) { return m.distance > 2; }), matches.end());
      std::vector<std::string> stringMatches;
      std::transform(matches.begin(), matches.end(), std::back_inserter(stringMatches), [](auto &m) { return m.name; });
      PRECICE_ERROR("The tag <{}> in the configuration contains an unknown attribute \"{}\". Did you mean \"{}\"?", _fullName, name, fmt::join(stringMatches, ", "));
    }
    PRECICE_ERROR("The tag <{}> in the configuration contains an unknown attribute \"{}\". Expected attributes are {}.", _fullName, name, fmt::join(expected, ", "));
  }

  // Read all attributes
  for (auto &attribute : _attributes) {
    std::visit(
        [&aAttributes](auto &attribute) { attribute.readValue(aAttributes); },
        attribute);
  }
}

std::vector<std::string> XMLTag::getAttributeNames() const
{
  std::vector<std::string> names;
  for (const auto &attribute : _attributes) {
    names.push_back(xml::getName(attribute));
  }
  return names;
}

void XMLTag::areAllSubtagsConfigured() const
{
  for (const auto &tag : _subtags) {
    std::string ns         = tag->_namespace;
    bool        configured = tag->isConfigured();

    bool occurOnce       = tag->getOccurrence() == OCCUR_ONCE;
    bool occurOnceOrMore = tag->getOccurrence() == OCCUR_ONCE_OR_MORE;

    if (not ns.empty()) {
      auto nsIter = _configuredNamespaces.find(ns);
      PRECICE_ASSERT(nsIter != _configuredNamespaces.end());
      configured |= nsIter->second;
    }

    if ((not configured) && (occurOnce || occurOnceOrMore)) {

      if (tag->getNamespace().empty()) {
        PRECICE_ERROR("Tag <{}> was not found but is required to occur at least once.", tag->getName());
      } else {
        PRECICE_ERROR("Tag <{}:... > was not found but is required to occur at least once.", tag->getNamespace());
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
  for (auto &attribute : _attributes) {
    std::visit(
        [](auto &attribute) { attribute.setRead(false); },
        attribute);
  }
  for (auto &tag : _subtags) {
    tag->_configured = false;
    tag->resetAttributes();
  }
}

std::string getName(const XMLTag::Attribute &attribute)
{
  return std::visit([](auto &attribute) { return attribute.getName(); }, attribute);
}

XMLTag getRootTag()
{
  static NoPListener listener;
  return XMLTag(listener, "configuration", XMLTag::OCCUR_ONCE);
}

void configure(
    XMLTag &                                  tag,
    const precice::xml::ConfigurationContext &context,
    std::string_view                          configurationFilename)
{
  logging::Logger _log("xml");
  PRECICE_TRACE(tag.getFullName(), configurationFilename);

  NoPListener nopListener;
  XMLTag      root(nopListener, "", XMLTag::OCCUR_ONCE);

  precice::xml::ConfigParser p(configurationFilename, context, std::make_shared<XMLTag>(tag));

  root.addSubtag(tag);
}

std::string_view XMLTag::getOccurrenceString(XMLTag::Occurrence occurrence)
{
  if (occurrence == XMLTag::OCCUR_ARBITRARY) {
    return "0..*";
  } else if (occurrence == XMLTag::OCCUR_NOT_OR_ONCE) {
    return "0..1";
  } else if (occurrence == XMLTag::OCCUR_ONCE) {
    return "1";
  } else if (occurrence == XMLTag::OCCUR_ONCE_OR_MORE) {
    return "1..*";
  }
  return "";
}
} // namespace precice::xml

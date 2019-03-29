#include "xml/XMLTag.hpp"
#include "utils/Helpers.hpp"
#include "utils/String.hpp"

namespace precice
{
namespace xml
{

XMLTag::XMLTag(
    Listener &         listener,
    const std::string &tagName,
    Occurrence         occurrence,
    const std::string &xmlNamespace)
    : _listener(listener),
      _name(tagName),
      _namespace(xmlNamespace),
      _occurrence(occurrence)
{
  if (not _namespace.empty()) {
    _fullName = _namespace + ":" + _name;
  } else {
    _fullName = _name;
  }
}

XMLTag& XMLTag::setDocumentation(const std::string &documentation)
{
  _doc = documentation;
  return *this;
}

XMLTag& XMLTag::addNamespace(const std::string &namespaceName)
{
  _namespaces.push_back(namespaceName);
  return *this;
}

XMLTag& XMLTag::addSubtag(const XMLTag &tag)
{
  TRACE(tag._fullName);
  assertion(tag._name != std::string(""));
  if (not tag._namespace.empty()) {
    _configuredNamespaces[tag._namespace] = false;
  }

  _subtags.push_back(std::make_shared<XMLTag>(tag));
  return *this;
}

XMLTag& XMLTag::addAttribute(const XMLAttribute<double> &attribute)
{
  TRACE(attribute.getName());
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _doubleAttributes.insert(std::pair<std::string, XMLAttribute<double>>(attribute.getName(), attribute));
  return *this;
}

XMLTag& XMLTag::addAttribute(const XMLAttribute<int> &attribute)
{
  TRACE(attribute.getName());
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _intAttributes.insert(std::pair<std::string, XMLAttribute<int>>(attribute.getName(), attribute));
  return *this;
}

XMLTag& XMLTag::addAttribute(const XMLAttribute<std::string> &attribute)
{
  TRACE(attribute.getName());
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _stringAttributes.insert(std::pair<std::string, XMLAttribute<std::string>>(attribute.getName(), attribute));
  return *this;
}

XMLTag& XMLTag::addAttribute(const XMLAttribute<bool> &attribute)
{
  TRACE(attribute.getName());
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _booleanAttributes.insert(std::pair<std::string, XMLAttribute<bool>>(attribute.getName(), attribute));
  return *this;
}

XMLTag& XMLTag::addAttribute(const XMLAttribute<Eigen::VectorXd> &attribute)
{
  TRACE(attribute.getName());
  assertion(not utils::contained(attribute.getName(), _attributes));
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
  assertion(iter != _doubleAttributes.end());
  return iter->second.getValue();
}

int XMLTag::getIntAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<int>>::const_iterator iter;
  iter = _intAttributes.find(name);
  assertion(iter != _intAttributes.end());
  return iter->second.getValue();
}

const std::string &XMLTag::getStringAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<std::string>>::const_iterator iter;
  iter = _stringAttributes.find(name);
  assertion(iter != _stringAttributes.end(), name);
  return iter->second.getValue();
}

bool XMLTag::getBooleanAttributeValue(const std::string &name) const
{
  std::map<std::string, XMLAttribute<bool>>::const_iterator iter;
  iter = _booleanAttributes.find(name);
  assertion(iter != _booleanAttributes.end());
  return iter->second.getValue();
}

Eigen::VectorXd XMLTag::getEigenVectorXdAttributeValue(const std::string &name, int dimensions) const
{
  TRACE(name, dimensions);
  // std::map<std::string, XMLAttribute<utils::DynVector> >::const_iterator iter;
  auto iter = _eigenVectorXdAttributes.find(name);
  assertion(iter != _eigenVectorXdAttributes.end());
  CHECK(iter->second.getValue().size() >= dimensions,
        "Vector attribute \"" << name << "\" of tag <" << getFullName()
                              << "> has less dimensions than required (" << iter->second.getValue().size()
                              << " instead of " << dimensions << ")!");

  // Read only first "dimensions" components of the parsed vector values
  Eigen::VectorXd        result(dimensions);
  const Eigen::VectorXd &parsed = iter->second.getValue();
  for (int i = 0; i < dimensions; i++) {
    result[i] = parsed[i];
  }
  DEBUG("Returning value = " << result);
  return result;
}

// new readattributes fx
void XMLTag::readAttributes(std::map<std::string, std::string> &aAttributes)
{
  TRACE();

  for (auto &element : aAttributes) {
    auto name = element.first;

    if (not utils::contained(name, _attributes)) {
      ERROR("Wrong attribute \"" << name << '\"');
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
  TRACE();
//  using utils::contained;
//  std::set<std::string> readNames;
  for (int i=0; i < xmlReader->getAttributeCount(); i++){
    std::string name = xmlReader->getAttributeName(i);
    if (not utils::contained(name, _attributes)){
      std::string error = "Wrong attribute \"" + name + "\"";
      throw error;
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
//      throw "Internal error in readAttributes";
//    }
//    readNames.insert(name);
  }

//  // Check if all attributes are read
//  for (const std::string& name : _attributes){
//    if (not contained(name, readNames)){
//
//      std::ostringstream stream;
//      stream << "Attribute \"" << name << "\" is not defined";
//      throw stream.str();
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
      assertion(utils::contained(ns, _configuredNamespaces));
      configured |= _configuredNamespaces.find(ns)->second;
    }

    if ((not configured) && (occurOnce || occurOnceOrMore)) {

      if (tag->getNamespace().empty()) {
        ERROR("Tag <" << tag->getName() << "> is missing");
      } else {
        ERROR("Tag <" << tag->getNamespace() << ":...> is missing");
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

  for (auto tag : _subtags) {
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

std::string XMLTag::printDTD(const bool start) const
{
  std::ostringstream dtd;

  if (start)
    dtd << "<!DOCTYPE " << _fullName << " [\n";

  dtd << "<!ELEMENT " << _fullName << " ";

  if (not _subtags.empty()) {

    dtd << "(";

    bool first = true;
    for (auto const subtag : _subtags) {

      std::string occurrenceChar = "";

      Occurrence occ = subtag->getOccurrence();

      if (occ == OCCUR_ARBITRARY)
        occurrenceChar = "*";
      else if (occ == OCCUR_NOT_OR_ONCE)
        occurrenceChar = "?";
      else if (occ == OCCUR_ONCE_OR_MORE)
        occurrenceChar = "+";

      dtd << (first ? "" : ", ") << subtag->getFullName() << occurrenceChar;
      first = false;
    }

    dtd << ")>\n";
  } else {
    dtd << "EMPTY>\n";
  }

  for (const auto &pair : _doubleAttributes) {
    dtd << pair.second.printDTD(_fullName);
  }

  for (const auto &pair : _intAttributes) {
    dtd << pair.second.printDTD(_fullName);
  }

  for (const auto &pair : _stringAttributes) {
    dtd << pair.second.printDTD(_fullName);
  }

  for (const auto &pair : _booleanAttributes) {
    dtd << pair.second.printDTD(_fullName);
  }

  for (const auto &pair : _eigenVectorXdAttributes) {
    dtd << pair.second.printDTD(_fullName);
  }

  if (not _subtags.empty()) {
    for (auto const subtag : _subtags) {
      dtd << subtag->printDTD();
    }
  }

  dtd << '\n';

  if (start)
    dtd << "]>\n";

  return dtd.str();
}

std::string XMLTag::printDocumentation(int indentation) const
{
  TRACE(indentation);
  const int   linewidth = 1000;
  std::string indent;
  for (int i = 0; i < indentation; i++) {
    indent += " ";
  }

  std::ostringstream doc;
  doc << indent << "<!-- TAG " << _fullName << '\n';
  if (not _doc.empty()) {
    std::string indentedDoc = indent + "         " + _doc;
    doc << utils::wrapText(indentedDoc, linewidth, indentation + 9);
    doc << '\n';
  }
  doc << indent << "         (can occur " << getOccurrenceString(_occurrence) << " times)";

  for (const auto &pair : _doubleAttributes) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : _intAttributes) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : _stringAttributes) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : _booleanAttributes) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : _eigenVectorXdAttributes) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  doc << " -->\n";
  std::ostringstream tagHead;
  tagHead << indent << "<" << _fullName;

  // Print XML namespaces, necessary for correct XML format and display in browser
  for (const std::string &namespaceName : _namespaces) {
    tagHead << " xmlns:" << namespaceName << "=\"precice." << namespaceName << "\"";
  }

  for (const auto &pair : _doubleAttributes) {
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  for (const auto &pair : _intAttributes) {
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  for (const auto &pair : _stringAttributes) {
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  for (const auto &pair : _booleanAttributes) {
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  for (const auto &pair : _eigenVectorXdAttributes) {
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  doc << utils::wrapText(tagHead.str(), linewidth, indentation + 3);

  if (not _subtags.empty()) {
    doc << ">\n\n";
    for (auto const subtag : _subtags) {
      doc << subtag->printDocumentation(indentation + 3);
    }
    doc << indent << "</" << _fullName << ">\n\n";
  } else {
    doc << "/>\n\n";
  }

  return doc.str();
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
    XMLTag &           tag,
    const std::string &configurationFilename)
{
  logging::Logger _log("xml");
  TRACE(tag.getFullName(), configurationFilename);

  NoPListener nopListener;
  XMLTag      root(nopListener, "", XMLTag::OCCUR_ONCE);

  precice::xml::ConfigParser p(configurationFilename, std::make_shared<XMLTag>(tag));

  root.addSubtag(tag);
}

std::string XMLTag::getOccurrenceString(Occurrence occurrence) const
{
  if (occurrence == OCCUR_ARBITRARY) {
    return std::string("0..*");
  } else if (occurrence == OCCUR_NOT_OR_ONCE) {
    return std::string("0..1");
  } else if (occurrence == OCCUR_ONCE) {
    return std::string("1");
  } else if (occurrence == OCCUR_ONCE_OR_MORE) {
    return std::string("1..*");
  }
  ERROR("Unknown occurrence type = " << occurrence);
  return "";
}
}
} // namespace precice, xml

//std::ostream& operator<<
//(
//  std::ostream&                 os,
//  const precice::xml::XMLTag& tag )
//{
//  os << tag.printDocumentation(80, 0);
//  return os;
//}

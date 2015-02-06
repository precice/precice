// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "utils/xml/XMLTag.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/String.hpp"

namespace precice {
namespace utils {

tarch::logging::Log precice::utils::XMLTag:: _log ("precice::utils::XMLTag");

XMLTag:: XMLTag
(
  Listener&          listener,
  const std::string& tagName,
  Occurrence         occurrence,
  const std::string& xmlNamespace )
:
  _listener(listener),
  _name(tagName),
  _namespace(xmlNamespace),
  _fullName(),
  _configured(false),
  _occurrence(occurrence),
  _namespaces(),
  _subtags(),
  _configuredNamespaces(),
  _attributes(),
  _doubleAttributes(),
  _intAttributes(),
  _stringAttributes(),
  _booleanAttributes(),
  _vector2DAttributes(),
  _vector3DAttributes(),
  _dynVectorAttributes()
{
  if (not _namespace.empty()){
    _fullName = _namespace + ":" + _name;
  }
  else {
    _fullName = _name;
  }
}

void XMLTag:: setDocumentation
(
  const std::string& documentation )
{
  _doc = documentation;
}

void XMLTag:: addNamespace
(
  const std::string& namespaceName )
{
  _namespaces.push_back(namespaceName);
}

void XMLTag:: addSubtag
(
  const XMLTag& tag )
{
  preciceTrace1("addSubtag()", tag._fullName);
  assertion (tag._name != std::string(""));
  if (not tag._namespace.empty()){
    _configuredNamespaces[tag._namespace] = false;
  }

  XMLTag* copy = new XMLTag(tag);   // resolves mingw problem
  _subtags.push_back(copy);
}

//void XMLTag:: removeSubtag
//(
//  const std::string& tagName )
//{
//  std::vector<XMLTag*>::iterator iter;
//  for ( iter=_subtags.begin(); iter != _subtags.end(); iter++ ){
//    if ((*iter)->getName() == tagName) {
//      delete *iter;  // MARK Bernhard, mingw
//      _subtags.erase(iter);
//      return;
//    }
//  }
//  preciceError ( "removeAttribute()", "Subtag \"" << tagName << "\" does not exist!" );
//}

void XMLTag:: addAttribute
(
  const XMLAttribute<double>& attribute )
{
  preciceTrace1 ( "addAttribute<double>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _doubleAttributes.insert(std::pair<std::string,XMLAttribute<double> >
                           (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<int>& attribute )
{
  preciceTrace1 ( "addAttribute<int>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _intAttributes.insert ( std::pair<std::string,
                          XMLAttribute<int> > (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<std::string>& attribute )
{
  preciceTrace1 ( "addAttribute<string>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _stringAttributes.insert ( std::pair<std::string,XMLAttribute<std::string> >
                             (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<bool>& attribute )
{
  preciceTrace1 ( "addAttribute<bool>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _booleanAttributes.insert ( std::pair<std::string,XMLAttribute<bool> >
      (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<utils::Vector2D>& attribute )
{
  preciceTrace1 ( "addAttribute<Vector2D>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _vector2DAttributes.insert ( std::pair<std::string,XMLAttribute<utils::Vector2D> >
      (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<utils::Vector3D>& attribute )
{
  preciceTrace1 ( "addAttribute<Vector3D>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _vector3DAttributes.insert(std::pair<std::string,XMLAttribute<utils::Vector3D> >
                             (attribute.getName(), attribute));
}

void XMLTag:: addAttribute
(
  const XMLAttribute<utils::DynVector>& attribute )
{
  preciceTrace1 ( "addAttribute<DynVector>()", attribute.getName() );
  assertion(not utils::contained(attribute.getName(), _attributes));
  _attributes.insert(attribute.getName());
  _dynVectorAttributes.insert (
      std::pair<std::string,XMLAttribute<utils::DynVector> >
      (attribute.getName(), attribute) );
}

bool XMLTag:: hasAttribute
(
  const std::string& attributeName )
{
//  bool has = false;
//  using utils::contained;
//  if (contained(attributeName, _intAttributes)) has = true;
//  else if (contained(attributeName, _doubleAttributes)) has = true;
//  else if (contained(attributeName, _stringAttributes)) has = true;
//  else if (contained(attributeName, _booleanAttributes)) has = true;
//  else if (contained(attributeName, _vector2DAttributes)) has = true;
//  else if (contained(attributeName, _vector3DAttributes)) has = true;
//  else if (contained(attributeName, _dynVectorAttributes)) has = true;
  return utils::contained(attributeName, _attributes);
}

//void XMLTag:: removeAttribute
//(
//  const std::string& attributeName )
//{
//  using utils::contained;
//  if (contained(attributeName, _intAttributes))
//    _intAttributes.erase (attributeName);
//  else if (contained(attributeName, _doubleAttributes))
//    _doubleAttributes.erase (attributeName);
//  else if (contained(attributeName, _stringAttributes))
//    _stringAttributes.erase (attributeName);
//  else if (contained(attributeName, _booleanAttributes))
//    _booleanAttributes.erase ( attributeName );
//  else if (contained(attributeName, _vector2DAttributes))
//    _vector2DAttributes.erase (attributeName);
//  else if (contained(attributeName, _vector3DAttributes))
//    _vector3DAttributes.erase (attributeName);
//  else if (contained(attributeName, _dynVectorAttributes))
//    _dynVectorAttributes.erase (attributeName);
//  else {
//    preciceError ( "removeAttribute()",
//                   "Attribute \"" << attributeName << "\" does not exist!" );
//  }
//}

//const XMLTag& XMLTag:: getTag
//(
//  const std::string& tagName ) const
//{
//  for (size_t i=0; i < _subtags.size(); i++) {
//    if (_subtags.at(i)->getName() == tagName)
//      return *_subtags.at(i);
//  }
//  preciceError ( "getTag(.)", "Tag with name " << tagName
//      << " does not exist for tag " << _name );
//}

void XMLTag:: parse
(
  XMLTag::XMLReader* xmlReader )
{
  preciceTrace1("parse()", _fullName);
  try {
    resetAttributes();
    if (xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT){
      assertion(xmlReader->getNodeName() != NULL);
      preciceDebug("reading attributes of tag " << xmlReader->getNodeName());
      readAttributes(xmlReader);
      _listener.xmlTagCallback(*this);
    }

    if (_subtags.size() > 0){
      while (xmlReader->read()){
        if (xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT){
          assertion(xmlReader->getNodeName() != NULL);
          preciceDebug("reading subtag " << xmlReader->getNodeName()
                       << " of tag " << _fullName);
          parseSubtag(xmlReader);
        }
        else if (xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT_END){
          assertion(xmlReader->getNodeName() != NULL);
          if (std::string(xmlReader->getNodeName()) == _fullName){
            preciceDebug("end of tag " << xmlReader->getNodeName());
            areAllSubtagsConfigured();
            _configured = true;
            _listener.xmlEndTagCallback(*this);
            //resetAttributes();
            return;
          }
          else {
            std::ostringstream stream;
            if (not _fullName.empty()){
              stream << "Invalid closing tag </" << xmlReader->getNodeName() << ">";
              throw stream.str();
            }
            //else {
            //  stream << "Invalid closing tag </" << xmlReader->getNodeName()
            //         << ">, expected closing of tag <" << _fullName << "> instead";
            //}
          }
        }
      }
      if (not _name.empty()){
        std::ostringstream error;
        error << "Missing closing tag </" << _fullName << ">";
        throw error.str();
      }
    }
    else {
      _configured = true;
    }
  }
  catch (std::string errorMsg){
    if (not _name.empty()){
      errorMsg += "\n   in tag <" + _fullName + ">";
    }
    throw errorMsg;
  }
}

void XMLTag:: parseSubtag
(
  XMLTag::XMLReader* xmlReader )
{
  preciceTrace1("parseSubtag()", _fullName);
  bool success = false;
  foreach (XMLTag* tag, _subtags){
    if (std::string(xmlReader->getNodeName()) == tag->getFullName()){
      if (tag->isConfigured()
          && (tag->getOccurrence() == OCCUR_ONCE
              || tag->getOccurrence() == OCCUR_NOT_OR_ONCE))
      {
        std::string error = "Tag <" + tag->getFullName() +
                            "> is not allowed to occur multiple times";
        throw error;
      }
      tag->parse(xmlReader);
      if (not tag->_namespace.empty()){
        assertion(utils::contained(tag->_namespace, _configuredNamespaces));
        _configuredNamespaces[tag->_namespace] = true;
      }
      success = true;
      break;   // since the xmlReader should not advance further
    }
  }
  if (not success){
    std::ostringstream stream;
    stream << "Wrong tag <" << xmlReader->getNodeName() << ">";
    throw stream.str();
  }

}

double XMLTag:: getDoubleAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<double> >::const_iterator iter;
  iter = _doubleAttributes.find(name);
  assertion(iter != _doubleAttributes.end());
  return iter->second.getValue();
}

int XMLTag:: getIntAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<int> >::const_iterator iter;
  iter = _intAttributes.find(name);
  assertion (iter  != _intAttributes.end());
  return iter->second.getValue();
}

const std::string& XMLTag:: getStringAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<std::string> >::const_iterator iter;
  iter = _stringAttributes.find(name);
  assertion1 (iter != _stringAttributes.end(), name);
  return iter->second.getValue();
}

bool XMLTag:: getBooleanAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<bool> >::const_iterator iter;
  iter = _booleanAttributes.find(name);
  assertion (iter  != _booleanAttributes.end());
  return iter->second.getValue();
}

const utils::Vector2D& XMLTag:: getVector2DAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<utils::Vector2D> >::const_iterator iter;
  iter = _vector2DAttributes.find(name);
  assertion (iter  != _vector2DAttributes.end());
  return iter->second.getValue();
}

const utils::Vector3D& XMLTag:: getVector3DAttributeValue
(
  const std::string& name ) const
{
  std::map<std::string,XMLAttribute<utils::Vector3D> >::const_iterator iter;
  iter = _vector3DAttributes.find(name);
  assertion (iter  != _vector3DAttributes.end());
  return iter->second.getValue();
}

utils::DynVector XMLTag:: getDynVectorAttributeValue
(
  const std::string& name,
  int                dimensions ) const
{
  preciceTrace2("getDynVectorAttributeValue()", name, dimensions);
  std::map<std::string,XMLAttribute<utils::DynVector> >::const_iterator iter;
  iter = _dynVectorAttributes.find(name);
  assertion (iter  != _dynVectorAttributes.end());
  preciceCheck(iter->second.getValue().size() >= dimensions,
               "getDynVectorAttributeValue()", "Vector attribute \""
               << name << "\" of tag <" << getFullName()
               << "> has less dimensions than required ("
               << iter->second.getValue().size() << " instead of "
               << dimensions << ")!" );

  // Read only first "dimensions" components of the parsed vector values
  utils::DynVector result(dimensions);
  const utils::DynVector& parsed = iter->second.getValue();
  for (int i=0; i < dimensions; i++){
    result[i] = parsed[i];
  }
  preciceDebug("Returning value = " << result);
  return result;
}

void XMLTag:: readAttributes
(
  XMLReader* xmlReader )
{
  preciceTrace("readAttributes()");
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
//  foreach (const std::string& name, _attributes){
//    if (not contained(name, readNames)){
//
//      std::ostringstream stream;
//      stream << "Attribute \"" << name << "\" is not defined";
//      throw stream.str();
//    }
//  }

  typedef std::map<std::string,XMLAttribute<double> >::value_type DoublePair;
  foreach (DoublePair & pair, _doubleAttributes){
     pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<int> >::value_type IntPair;
  foreach ( IntPair & pair, _intAttributes){
    pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<std::string> >::value_type StringPair;
  foreach ( StringPair & pair, _stringAttributes ){
    pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<bool> >::value_type BoolPair;
  foreach ( BoolPair & pair, _booleanAttributes ){
    pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector2D> >::value_type Vec2DPair;
  foreach ( Vec2DPair& pair, _vector2DAttributes ){
    pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector3D> >::value_type Vec3DPair;
  foreach ( Vec3DPair& pair, _vector3DAttributes ){
    pair.second.readValue(xmlReader);
  }

  typedef std::map<std::string,XMLAttribute<utils::DynVector> >::value_type DynVecPair;
  foreach ( DynVecPair& pair, _dynVectorAttributes ){
    pair.second.readValue(xmlReader);
  }
}

void XMLTag:: areAllSubtagsConfigured() const
{
  foreach (XMLTag* tag, _subtags){
    std::string ns = tag->_namespace;
    bool configured = tag->isConfigured();
    if (not ns.empty()){
      assertion(utils::contained(ns, _configuredNamespaces));
      configured |= _configuredNamespaces.find(ns)->second;
    }
    bool occurOnce = tag->getOccurrence() == OCCUR_ONCE;
    bool occurOnceOrMore = tag->getOccurrence() == OCCUR_ONCE_OR_MORE;

    if ( (not configured) && (occurOnce || occurOnceOrMore) ){
      std::ostringstream stream;
      if (tag->getNamespace().empty()){
        stream << "Tag <" << tag->getName() << "> is missing";
      }
      else {
        stream << "Tag <" << tag->getNamespace() << ":...> is missing";
      }
      throw stream.str();
    }
  }
}

void XMLTag:: resetAttributes()
{
  _configured = false;

  typedef std::map<std::string,bool>::value_type Pair;
  foreach (Pair& pair, _configuredNamespaces){
    pair.second = false;
  }

  typedef std::map<std::string,XMLAttribute<double> >::value_type DoublePair;
  foreach (DoublePair& pair, _doubleAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<int> >::value_type IntPair;
  foreach (IntPair& pair, _intAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<std::string> >::value_type StringPair;
  foreach (StringPair& pair, _stringAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<bool> >::value_type BoolPair;
  foreach (BoolPair& pair, _booleanAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector2D> >::value_type Vec2DPair;
  foreach (Vec2DPair& pair, _vector2DAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector3D> >::value_type Vec3DPair;
  foreach (Vec3DPair& pair, _vector3DAttributes){
    pair.second.setRead(false);
  }

  typedef std::map<std::string,XMLAttribute<utils::DynVector> >::value_type DynVecPair;
  foreach (DynVecPair& pair, _dynVectorAttributes){
    pair.second.setRead(false);
  }

  std::vector<XMLTag*>::iterator subtagIter;
  foreach (XMLTag* tag, _subtags){
    tag->_configured = false;
    tag->resetAttributes();
  }
}

void XMLTag:: clear()
{
   _doubleAttributes.clear();
   _intAttributes.clear();
   _stringAttributes.clear();
   _booleanAttributes.clear();
   _vector2DAttributes.clear();
   _vector3DAttributes.clear();
   _dynVectorAttributes.clear();
   _subtags.clear();
}

std::string XMLTag:: printDocumentation
(
  int linewidth,
  int indentation ) const
{
  preciceTrace1("printDocumentation()", indentation);
  std::string indent;
  for (int i=0; i < indentation; i++){
    indent += " ";
  }

  std::ostringstream doc;
  //std::ostringstream tagDoc;
  doc << indent << "<!-- TAG " << _fullName << std::endl;
  if (not _doc.empty()){
    std::string indentedDoc = indent + "         " + _doc;
    doc <<  wrapText(indentedDoc, linewidth, indentation + 9);
    doc << std::endl;
  }
  doc << indent << "         (can occur " << getOccurrenceString(_occurrence) << " times)";

//  if (not _attributes.empty()){
//    doc << std::endl << indent << "     Attribute documentation:";
//  }

  typedef std::map<std::string,XMLAttribute<double> >::value_type DoublePair;
  foreach (const DoublePair& pair, _doubleAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<int> >::value_type IntPair;
  foreach (const IntPair& pair, _intAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<std::string> >::value_type StringPair;
  foreach (const StringPair& pair, _stringAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<bool> >::value_type BoolPair;
  foreach (const BoolPair& pair, _booleanAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector2D> >::value_type Vec2DPair;
  foreach (const Vec2DPair& pair, _vector2DAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector3D> >::value_type Vec3DPair;
  foreach (const Vec3DPair& pair, _vector3DAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  typedef std::map<std::string,XMLAttribute<utils::DynVector> >::value_type DynVecPair;
  foreach (const DynVecPair& pair, _dynVectorAttributes){
    std::ostringstream attrDoc;
    doc << std::endl;
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  doc << " -->" << std::endl;
  std::ostringstream tagHead;
  tagHead << indent << "<" << _fullName;

  // Print XML namespaces, necessary for correct XML format and display in browser
  foreach (const std::string& namespaceName, _namespaces){
    tagHead << " xmlns:" << namespaceName << "=\"precice." << namespaceName << "\"";
  }

  typedef std::map<std::string,XMLAttribute<double> >::value_type DoublePair;
  foreach (const DoublePair& pair, _doubleAttributes){
    //std::ostringstream attrDoc;
    //doc << std::endl;
    tagHead << indent << "   " << pair.second.printDocumentation();
    //doc << wrapText(attrDoc.str(), linewidth, indentation + 6);
  }

  typedef std::map<std::string,XMLAttribute<int> >::value_type IntPair;
  foreach (const IntPair& pair, _intAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  typedef std::map<std::string,XMLAttribute<std::string> >::value_type StringPair;
  foreach (const StringPair& pair, _stringAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  typedef std::map<std::string,XMLAttribute<bool> >::value_type BoolPair;
  foreach (const BoolPair& pair, _booleanAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector2D> >::value_type Vec2DPair;
  foreach (const Vec2DPair& pair, _vector2DAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  typedef std::map<std::string,XMLAttribute<utils::Vector3D> >::value_type Vec3DPair;
  foreach (const Vec3DPair& pair, _vector3DAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  typedef std::map<std::string,XMLAttribute<utils::DynVector> >::value_type DynVecPair;
  foreach (const DynVecPair& pair, _dynVectorAttributes){
    tagHead << indent << "   " << pair.second.printDocumentation();
  }

  doc << wrapText(tagHead.str(), linewidth, indentation + 3);

  if (not _subtags.empty()){
    doc << ">" << std::endl << std::endl;
    foreach (const XMLTag* subtag, _subtags){
      doc << subtag->printDocumentation(linewidth, indentation + 3);
    }
    doc << indent << "</" << _fullName << ">" << std::endl << std::endl;
  }
  else {
    doc << "/>" << std::endl << std::endl;
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

void configure
(
  XMLTag&            tag,
  const std::string& configurationFilename )
{
  tarch::logging::Log _log("precice::utils");
  preciceTrace2("configure()", tag.getFullName(), configurationFilename);
  //bool success = false;
  tarch::irr::io::IrrXMLReader* xmlReader =
    tarch::irr::io::createIrrXMLReader(configurationFilename.c_str());
  preciceCheck(xmlReader != NULL, "configure()",
               "Could not create XML reader for file \"" << configurationFilename
               << "\"!");
  preciceCheck(xmlReader->read(), "configure()",
               "XML reader doesn't recognize a valid XML tag in file \""
               << configurationFilename << "\"!" );
  preciceCheck(xmlReader->getNodeType() != tarch::irr::io::EXN_NONE, "configure()",
               "XML reader found only invalid XML tag in file \""
               << configurationFilename << "\"!" );
  //bool foundTag = false;
  NoPListener nopListener;
  XMLTag root(nopListener, "", XMLTag::OCCUR_ONCE);
  root.addSubtag(tag);
  try {
    root.parse(xmlReader);
  }
  catch (std::string errorMsg){
    preciceError("configure()", "Parsing XML file \""
                 << configurationFilename << "\"" << std::endl << errorMsg);
  }

//  while (xmlReader->read()){
//    if (xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT){
//      if (root.getFullName() == xmlReader->getNodeName()){
//        foundTag = true;
//        try {
//          root.parse(xmlReader);
//          success = true;
//        }
//        catch (std::string errorMsg){
//          preciceError("configure()", "Parsing XML file \""
//                       << configurationFilename << "\"" << std::endl << errorMsg);
//        }
//      }
//      else {
//        preciceError("configure()", "Found wrong tag <" << xmlReader->getNodeName()
//                     << "> in XML file \"" << configurationFilename << "\"!");
//      }
//    }
//    else if (xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT_END){
//      preciceError("configure()", "Found wrong end tag </" << xmlReader->getNodeName()
//                           << "> in XML file \"" << configurationFilename << "\"!");
//    }
//  }
//  preciceCheck(foundTag, "configure()", "Did not find root tag <" << root.getFullName()
//               << "> in XML configuration \"" << configurationFilename << "\"!");
//  return success;
}

std::string XMLTag:: getOccurrenceString ( Occurrence occurrence ) const
{
  if (occurrence == OCCUR_ARBITRARY){
    return std::string("0..*");
  }
  else if (occurrence == OCCUR_NOT_OR_ONCE){
    return std::string("0..1");
  }
  else if (occurrence == OCCUR_ONCE){
    return std::string("1");
  }
  else if (occurrence == OCCUR_ONCE_OR_MORE){
    return std::string("1..*");
  }
  preciceError("getOccurrenceString()", "Unknown occurrence type = " << occurrence);
  return "";
}

}} // namespace precice, utils

//std::ostream& operator<<
//(
//  std::ostream&                 os,
//  const precice::utils::XMLTag& tag )
//{
//  os << tag.printDocumentation(80, 0);
//  return os;
//}

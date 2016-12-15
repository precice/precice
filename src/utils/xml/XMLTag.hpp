#pragma once

#include "XMLAttribute.hpp"
#include "utils/Globals.hpp"
#include "tarch/irr/XML.h"
#include <string>
#include <map>
#include <vector>
#include <set>
#include "logging/Logger.hpp"

namespace precice {
namespace utils {


/// Represents an XML tag to be configured automatically.
class XMLTag
{
public:

  typedef tarch::irr::io::IrrXMLReader XMLReader;

  /// Callback interface for configuration classes using XMLTag.
  struct Listener
  {
    virtual ~Listener() {};
    /**
     * @brief Callback at begin of XML tag.
     *
     * At this callback, the attributes of the callingTag are already parsed and
     * available, while the subtags are not yet parsed.
     */
    virtual void xmlTagCallback ( XMLTag& callingTag ) =0;

    /**
     * @brief Callback at end of XML tag and at end of subtag.
     *
     * At this callback, the attributes and all subtags of callingTag are parsed.
     * This callback is first done for the listener, and then for the parent tag
     * listener (if existing).
     */
    virtual void xmlEndTagCallback ( XMLTag& callingTag ) =0;
  };

  /// Types of occurrences of an XML tag.
  enum Occurrence {
    OCCUR_NOT_OR_ONCE,
    OCCUR_ONCE,
    OCCUR_ONCE_OR_MORE,
    OCCUR_ARBITRARY,
    OCCUR_ARBITRARY_NESTED
  };

  /**
   * @brief Standard constructor
   *
   * @param[in] listener Configuration listener receiving callbacks from this tag.
   * @param[in] name Name of the XML tag.
   * @param[in] occurrence Defines the occurrences of the tag in the configuration.
   * @param[in] xmlNamespace Defines a prefix/namespace for the tag. Tags with equal namespace or treated as group.
   */
  XMLTag (
    Listener&          listener,
    const std::string& name,
    Occurrence         Occurrence,
    const std::string& xmlNamespace = "" );

  /**
   * @brief Adds a description of the purpose of this XML tag.
   *
   * The description and more information is printed with printDocumentation().
   */
  void setDocumentation(const std::string& documentation);

  /**
   * @brief Adds a namespace to the tag.
   *
   * Only used for outputting correct XML format, such that, e.g., internet
   * browsers display no errors when viewing an XML configuration.
   */
  void addNamespace( const std::string& namespaceName );

  /// Adds an XML tag as subtag by making a copy of the given tag.
  void addSubtag ( const XMLTag& tag );

  /// Removes the XML subtag with given name
  //void removeSubtag ( const std::string& tagName );

  /// Adds a XML attribute by making a copy of the given attribute.
  void addAttribute ( const XMLAttribute<double>& attribute );

  // Adds a XML attribute by making a copy of the given attribute.
  void addAttribute ( const XMLAttribute<int>& attribute );

  /// Adds a XML attribute by making a copy of the given attribute.
  void addAttribute ( const XMLAttribute<std::string>& attribute );

  /// Adds a XML attribute by making a copy of the given attribute.
  void addAttribute ( const XMLAttribute<bool>& attribute );

  /// Adds a XML attribute by making a copy of the given attribute.
  void addAttribute ( const XMLAttribute<Eigen::VectorXd>& attribute );
  
  bool hasAttribute ( const std::string& attributeName );

  //void removeAttribute ( const std::string& attributeName );

  /**
   * @brief Returns name (without namespace).
   *
   * The method getFullName() returns the name with namespace.
   */
  const std::string& getName() const { return _name; }

  /// Returns xml namespace.
  const std::string& getNamespace() const { return _namespace; }

  /**
   * @brief Returns full name consisting of xml namespace + ":" + name.
   */
  const std::string& getFullName() const { return _fullName; }

  //const XMLTag& getTag ( const std::string& tagName ) const;

  double getDoubleAttributeValue ( const std::string& name ) const;

  int getIntAttributeValue ( const std::string& name ) const;

  const std::string& getStringAttributeValue ( const std::string& name ) const;

  bool getBooleanAttributeValue ( const std::string& name ) const;

  /**
   * @brief Returns Eigen vector attribute value with given dimensions.
   *
   * If the parsed vector has less dimensions then required, an error message
   * is thrown.
   *
   * @param[in] name Name of attribute.
   * @param[in] dimensions Dimensions of the vector to be returned.
   */
  Eigen::VectorXd getEigenVectorXdAttributeValue (
    const std::string& name,
    int                dimensions ) const;
  
  /// Parses the information from the xmlReader and calls XMLListener.
  void parse ( XMLReader* xmlReader );

  bool isConfigured() const { return _configured; }

  Occurrence getOccurrence() const { return _occurrence; }

  /// Removes all attributes and subtags
  void clear();

  /// Prints a documentation string for this tag.
  std::string printDocumentation (
    int linewidth,
    int indentation ) const;

private:

  static logging::Logger _log;

  Listener& _listener;

  /// Name of the tag.
  std::string _name;

  /// XML namespace of the tag.
  std::string _namespace;

  /// Combination of name and namespace: _namespace + ":" + _name
  std::string _fullName;

  std::string _doc;

  bool _configured;

  Occurrence _occurrence;

  std::vector<std::string> _namespaces;

  std::vector<XMLTag*> _subtags;

  std::map<std::string,bool> _configuredNamespaces;

  std::set<std::string> _attributes;

  std::map<std::string,XMLAttribute<double> > _doubleAttributes;

  std::map<std::string,XMLAttribute<int> > _intAttributes;

  std::map<std::string,XMLAttribute<std::string> > _stringAttributes;

  std::map<std::string,XMLAttribute<bool> > _booleanAttributes;

  std::map<std::string,XMLAttribute<Eigen::VectorXd> > _eigenVectorXdAttributes;

  void readAttributes ( XMLReader* xmlReader );

  void parseSubtag ( XMLReader*  xmlReader );

  void areAllSubtagsConfigured() const;

  void resetAttributes();

  std::string getOccurrenceString ( Occurrence occurrence ) const;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

/// No operation listener for tests.
struct NoPListener : public XMLTag::Listener
{
  virtual void xmlTagCallback ( XMLTag& callingTag ) {}
  virtual void xmlEndTagCallback ( XMLTag& callingTag ) {}
};

/**
 * @brief Returns an XMLTag::Listener that does nothing on callbacks.
 *
 * This is useful for tests, when the root tag to be specified in
 */
//NoPListener& getNoPListener();

/**
 * @brief Returns an empty root tag with name "configuration".
 *
 * A static NoPListener is added, and the occurrence is set to OCCUR_ONCE.
 */
XMLTag getRootTag();

/**
 * @brief Configures the given configuration from file configurationFilename.
 */
void configure (
  XMLTag&            tag,
  const std::string& configurationFilename );

///**
// * @brief Configures the given configuration from file configurationFilename.
// */
//template< typename CONFIG_T >
//bool configure
//(
//  CONFIG_T&          configuration,
//  const std::string& configurationFilename )
//{
//  logging::Logger _log ( "precice::utils" );
//  bool success = false;
//  tarch::irr::io::IrrXMLReader * xmlReader =
//    tarch::irr::io::createIrrXMLReader ( configurationFilename.c_str() );
//  preciceCheck ( xmlReader != NULL, "configure()",
//                 "Could not create XML reader for file \"" << configurationFilename
//                 << "\"!" );
//  preciceCheck ( xmlReader->read(), "configure()",
//                 "XML reader doesn't recognize a valid XML tag in file \""
//                 << configurationFilename << "\"!" );
//  preciceCheck ( xmlReader->getNodeType() != tarch::irr::io::EXN_NONE, "configure()",
//                 "XML reader found only invalid XML tag in file \""
//                 << configurationFilename << "\"!" );
//  bool foundTag = false;
//  while( xmlReader->read() ) {
//    if ( xmlReader->getNodeType() == tarch::irr::io::EXN_ELEMENT ) {
//      if ( CONFIG_T::getTag() == xmlReader->getNodeName() ) {
//        foundTag = true;
//        success |= configuration.parseSubtag ( xmlReader );
//      }
//    }
//  }
//  preciceCheck ( foundTag, "configure()", "Did not find suitable XML tag!" );
//  return success;
//};
//
///**
// * @brief Configures the given confugration from file configurationFilename.
// */
//template< typename CONFIG_T >
//bool configure
//(
//  std::shared_ptr<CONFIG_T> configuration,
//  const std::string&          configurationFilename )
//{
//  assertion ( configuration.get() != NULL );
//  return configure ( *configuration, configurationFilename );
//}

}} // namespace precice, utils

/**
 * @brief Adds documentation of tag to output stream os.
 */
//std::ostream& operator<< (
//  std::ostream&                 os,
//  const precice::utils::XMLTag& tag );


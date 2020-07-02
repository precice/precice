#pragma once

#include <Eigen/Core>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace xml {
class ConfigParser;
}
} // namespace precice

namespace precice {
namespace xml {

/// Tightly coupled to the parameters of SolverInterface()
struct ConfigurationContext {
  std::string name;
  int         rank;
  int         size;
};

/// Represents an XML tag to be configured automatically.
class XMLTag {
  friend class precice::xml::ConfigParser;

public:
  using Namespaces = typename std::vector<std::string>;

  using Subtags = typename std::vector<std::shared_ptr<XMLTag>>;

  template <typename T>
  using AttributeMap = typename std::map<std::string, XMLAttribute<T>>;

  /// Callback interface for configuration classes using XMLTag.
  struct Listener {

    Listener &operator=(Listener &&) = delete;

    virtual ~Listener(){};
    /**
     * @brief Callback at begin of XML tag.
     *
     * At this callback, the attributes of the callingTag are already parsed and
     * available, while the subtags are not yet parsed.
     */
    virtual void xmlTagCallback(ConfigurationContext const &context, XMLTag &callingTag) = 0;

    /**
     * @brief Callback at end of XML tag and at end of subtag.
     *
     * At this callback, the attributes and all subtags of callingTag are parsed.
     * This callback is first done for the listener, and then for the parent tag
     * listener (if existing).
     */
    virtual void xmlEndTagCallback(ConfigurationContext const &context, XMLTag &callingTag) = 0;
  };

  /// Types of occurrences of an XML tag.
  enum Occurrence {
    OCCUR_NOT_OR_ONCE,
    OCCUR_ONCE,
    OCCUR_ONCE_OR_MORE,
    OCCUR_ARBITRARY,
    OCCUR_ARBITRARY_NESTED
  };

  static std::string getOccurrenceString(Occurrence occurrence);

  /**
   * @brief Standard constructor
   *
   * @param[in] listener Configuration listener receiving callbacks from this tag.
   * @param[in] name Name of the XML tag.
   * @param[in] occurrence Defines the occurrences of the tag in the configuration.
   * @param[in] xmlNamespace Defines a prefix/namespace for the tag. Tags with equal namespace or treated as group.
   */
  XMLTag(
      Listener &  listener,
      std::string name,
      Occurrence  occurrence,
      std::string xmlNamespace = "");

  /**
   * @brief Adds a description of the purpose of this XML tag.
   *
   * The description and more information is printed with printDocumentation().
   */
  XMLTag &setDocumentation(const std::string &documentation);

  std::string getDocumentation() const
  {
    return _doc;
  };

  /**
   * @brief Adds a namespace to the tag.
   *
   * Only used for outputting correct XML format, such that, e.g., internet
   * browsers display no errors when viewing an XML configuration.
   */
  XMLTag &addNamespace(const std::string &namespaceName);

  /// Adds an XML tag as subtag by making a copy of the given tag.
  XMLTag &addSubtag(const XMLTag &tag);

  const Subtags &getSubtags() const
  {
    return _subtags;
  };

  /// Removes the XML subtag with given name
  //XMLTag& removeSubtag ( const std::string& tagName );

  /// Adds a XML attribute by making a copy of the given attribute.
  XMLTag &addAttribute(const XMLAttribute<double> &attribute);

  // Adds a XML attribute by making a copy of the given attribute.
  XMLTag &addAttribute(const XMLAttribute<int> &attribute);

  /// Adds a XML attribute by making a copy of the given attribute.
  XMLTag &addAttribute(const XMLAttribute<std::string> &attribute);

  /// Adds a XML attribute by making a copy of the given attribute.
  XMLTag &addAttribute(const XMLAttribute<bool> &attribute);

  /// Adds a XML attribute by making a copy of the given attribute.
  XMLTag &addAttribute(const XMLAttribute<Eigen::VectorXd> &attribute);

  bool hasAttribute(const std::string &attributeName);

  /**
   * @brief Returns name (without namespace).
   *
   * The method getFullName() returns the name with namespace.
   */
  const std::string &getName() const
  {
    return _name;
  }

  /// Returns xml namespace.
  const std::string &getNamespace() const
  {
    return _namespace;
  }

  const Namespaces &getNamespaces() const
  {
    return _namespaces;
  };

  /// Returns full name consisting of xml namespace + ":" + name.
  const std::string &getFullName() const
  {
    return _fullName;
  }

  double getDoubleAttributeValue(const std::string &name) const;

  int getIntAttributeValue(const std::string &name) const;

  const std::string &getStringAttributeValue(const std::string &name) const;

  bool getBooleanAttributeValue(const std::string &name) const;

  const AttributeMap<double> &getDoubleAttributes() const
  {
    return _doubleAttributes;
  };

  const AttributeMap<int> &getIntAttributes() const
  {
    return _intAttributes;
  };

  const AttributeMap<std::string> &getStringAttributes() const
  {
    return _stringAttributes;
  };

  const AttributeMap<bool> &getBooleanAttributes() const
  {
    return _booleanAttributes;
  };

  const AttributeMap<Eigen::VectorXd> &getEigenVectorXdAttributes() const
  {
    return _eigenVectorXdAttributes;
  };

  /**
   * @brief Returns Eigen vector attribute value with given dimensions.
   *
   * If the parsed vector has less dimensions then required, an error message
   * is thrown.
   *
   * @param[in] name Name of attribute.
   * @param[in] dimensions Dimensions of the vector to be returned.
   */
  Eigen::VectorXd getEigenVectorXdAttributeValue(
      const std::string &name,
      int                dimensions) const;

  bool isConfigured() const
  {
    return _configured;
  }

  Occurrence getOccurrence() const
  {
    return _occurrence;
  }

  /// Removes all attributes and subtags
  void clear();

  /// reads all attributes of this tag
  void readAttributes(const std::map<std::string, std::string> &aAttributes);

private:
  mutable logging::Logger _log{"xml::XMLTag"};

  Listener &_listener;

  /// Name of the tag.
  std::string _name;

  /// XML namespace of the tag.
  std::string _namespace;

  /// Combination of name and namespace: _namespace + ":" + _name
  std::string _fullName;

  std::string _doc;

  bool _configured = false;

  Occurrence _occurrence;

  Namespaces _namespaces;

  Subtags _subtags;

  std::map<std::string, bool> _configuredNamespaces;

  std::set<std::string> _attributes;

  AttributeMap<double> _doubleAttributes;

  AttributeMap<int> _intAttributes;

  AttributeMap<std::string> _stringAttributes;

  AttributeMap<bool> _booleanAttributes;

  AttributeMap<Eigen::VectorXd> _eigenVectorXdAttributes;

  void areAllSubtagsConfigured() const;

  void resetAttributes();
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

/// No operation listener for tests.
struct NoPListener : public XMLTag::Listener {
  void xmlTagCallback(ConfigurationContext const &context, XMLTag &callingTag) override {}
  void xmlEndTagCallback(ConfigurationContext const &context, XMLTag &callingTag) override {}
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

/// Configures the given configuration from file configurationFilename.
void configure(
    XMLTag &                                  tag,
    const precice::xml::ConfigurationContext &context,
    const std::string &                       configurationFilename);

} // namespace xml
} // namespace precice

/**
 * @brief Adds documentation of tag to output stream os.
 */
//std::ostream& operator<< (
//  std::ostream&                 os,
//  const precice::xml::XMLTag& tag );

#pragma once

#include <Eigen/Core>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include "logging/Logger.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::xml {
class ConfigParser;
}

namespace precice::xml {

/// Tightly coupled to the parameters of Participant()
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

  using Attribute = std::variant<
      XMLAttribute<double>,
      XMLAttribute<int>,
      XMLAttribute<std::string>,
      XMLAttribute<bool>,
      XMLAttribute<Eigen::VectorXd>>;

  using Attributes = typename std::vector<Attribute>;

  /// Callback interface for configuration classes using XMLTag.
  struct Listener {

    Listener &operator=(Listener &&) = delete;

    virtual ~Listener() = default;
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

  static std::string_view getOccurrenceString(Occurrence occurrence);

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
  XMLTag &setDocumentation(std::string_view documentation);

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

  /// Adds a hint for missing attributes, which will be displayed along the error message.
  void addAttributeHint(std::string name, std::string message);

  bool hasAttribute(const std::string &attributeName) const;

  template <typename Container>
  void addSubtags(const Container &subtags)
  {
    std::for_each(subtags.begin(), subtags.end(), [this](auto &s) { this->addSubtag(s); });
  }

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

  double getDoubleAttributeValue(const std::string &name, std::optional<double> default_value = std::nullopt) const;

  int getIntAttributeValue(const std::string &name, std::optional<int> default_value = std::nullopt) const;

  std::string getStringAttributeValue(const std::string &name, std::optional<std::string> default_value = std::nullopt) const;

  bool getBooleanAttributeValue(const std::string &name, std::optional<bool> default_value = std::nullopt) const;

  Eigen::VectorXd getEigenVectorXdAttributeValue(const std::string &name) const;

  std::vector<std::string> getAttributeNames() const;

  const Attributes &getAttributes() const
  {
    return _attributes;
  };

  bool isConfigured() const
  {
    return _configured;
  }

  Occurrence getOccurrence() const
  {
    return _occurrence;
  }

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

  Attributes _attributes;

  std::map<std::string, std::string> _attributeHints;

  void areAllSubtagsConfigured() const;

  void resetAttributes();
};

/// Returns the name of an Attribute
std::string getName(const XMLTag::Attribute &attribute);

// ------------------------------------------------------ HEADER IMPLEMENTATION

/// No operation listener for tests.
struct NoPListener : public XMLTag::Listener {
  void xmlTagCallback(ConfigurationContext const &context, XMLTag &callingTag) override {}
  void xmlEndTagCallback(ConfigurationContext const &context, XMLTag &callingTag) override {}
};

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
    std::string_view                          configurationFilename);

} // namespace precice::xml

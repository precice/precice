#pragma once

#include <libxml/SAX.h>
#include <map>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice
{
namespace xml
{
class XMLTag; // forward declaration to resolve circular import

class ConfigParser
{
public:
  /// Struct holding the read tag from xml file
  struct CTag {
    CTag()
        : m_Used(false) {}

    std::string m_Name;
    std::string m_Prefix;
    bool        m_Used;

    typedef std::map<std::string, std::string> AttributePair;
    AttributePair                              m_aAttributes;
    std::vector<CTag *>                        m_aSubTags;
  };

private:
  static precice::logging::Logger _log;

  std::vector<CTag *> m_AllTags;
  std::vector<CTag *> m_CurrentTags;

  precice::xml::XMLTag *m_pXmlTag;

  static void GenericErrorFunc(void *ctx, const char *msg, ...);

  /// Opens file and starts parsing
  void init(const std::string &filePath);

public:
  typedef CTag::AttributePair AttributePair;

  /// Parser ctor for Callback init
  ConfigParser(const std::string &filePath, XMLTag *pXmlTag);

  /// Parser ctor without Callbacks
  ConfigParser(const std::string &filePath);

  /// Removes all used tags
  ~ConfigParser();

  /// Reads the xml file
  int readXmlFile(FILE *f);

  /// Creates the handler with callbacks for the SAX interface
  xmlSAXHandler makeSaxHandler();

  /// Returns the root tag
  CTag *getRootTag();

  /**
   * @brief Connects the actual tags of an xml layer with the predefined tags
   * @param DefTags predefined tags
   * @param SubTags actual tags from xml file
   */
  void connectTags(std::vector<precice::xml::XMLTag *> &DefTags, std::vector<CTag *> &SubTags);

  /// Callback for Start-Tag
  static void OnStartElementNs(
      void *          ctx,
      const xmlChar * localname,
      const xmlChar * prefix,
      const xmlChar * URI,
      int             nb_namespaces,
      const xmlChar **namespaces,
      int             nb_attributes,
      int             nb_defaulted,
      const xmlChar **attributes);

  /// Callback for End-Tag
  static void OnEndElementNs(
      void *         ctx,
      const xmlChar *localname,
      const xmlChar *prefix,
      const xmlChar *URI);

  // Callback for text sections in xml file
  static void OnCharacters(void *ctx, const xmlChar *ch, int len);
};
}
}

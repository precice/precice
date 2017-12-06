#pragma once

#include <iostream>
#include <libxml/SAX.h>
#include <map>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "utils/xml/XMLTag.hpp"

namespace precice
{
namespace utils
{
class XMLTag;
}
}

namespace precice
{
namespace xml
{

class Parser
{
public:
  /// struct holding the read tag from xml file
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

  precice::utils::XMLTag *m_pXmlTag;

  static void GenericErrorFunc(void *ctx, const char *msg, ...);

  /// opens file and starts parsing
  void init(const std::string &filePath);

public:
  typedef CTag::AttributePair AttributePair;

  /// Parser ctor for Callback init
  Parser(const std::string &filePath, precice::utils::XMLTag *pXmlTag);

  /// Parser ctor without Callbacks
  Parser(const std::string &filePath);

  /// removes all used tags
  ~Parser();

  /// reads the xml file
  int readXmlFile(FILE *f);

  /// creates the handler with callbacks for the SAX interface
  xmlSAXHandler makeSaxHandler();

  /// returns the root tag
  CTag *getRootTag();

  /**
	 * @brief Connects the actual tags of an xml layer with the predefined tags
	 * @param DefTags predefined tags
	 * @param SubTags actual tags from xml file
	 */
  void connectTags(std::vector<precice::utils::XMLTag *> &DefTags, std::vector<CTag *> &SubTags);

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

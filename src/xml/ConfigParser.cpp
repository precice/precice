#include "xml/ConfigParser.hpp"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iterator>
#include <libxml/SAX.h>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace xml {

std::string decodeXML(std::string xml)
{
  static const std::map<std::string, char> escapes{{"&lt;", '<'}, {"&gt;", '>'}, {"&amp;", '&'}, {"&quot;", '"'}, {"&apos;", '\''}};
  while (true) {
    bool changes{false};
    for (const auto &kv : escapes) {
      auto position = xml.find(kv.first);
      if (position != std::string::npos) {
        xml.replace(position, kv.first.length(), 1, kv.second);
        changes = true;
      }
    }
    if (!changes) {
      break;
    }
  };
  return xml;
}

// ------------------------- Callback functions for libxml2  -------------------------

void OnStartElementNs(
    void *          ctx,
    const xmlChar * localname,
    const xmlChar * prefix,
    const xmlChar * URI,
    int             nb_namespaces,
    const xmlChar **namespaces,
    int             nb_attributes,
    int             nb_defaulted,
    const xmlChar **attributes)
{
  ConfigParser::CTag::AttributePair attributesMap;
  unsigned int                      index = 0;
  for (int indexAttribute = 0; indexAttribute < nb_attributes; ++indexAttribute, index += 5) {
    std::string attributeName(reinterpret_cast<const char *>(attributes[index]));

    auto        valueBegin = reinterpret_cast<const char *>(attributes[index + 3]);
    auto        valueEnd   = reinterpret_cast<const char *>(attributes[index + 4]);
    std::string value(valueBegin, valueEnd);

    attributesMap[attributeName] = decodeXML(value);
  }

  auto pParser = static_cast<ConfigParser *>(ctx);

  std::string sPrefix(prefix == nullptr ? "" : reinterpret_cast<const char *>(prefix));

  pParser->OnStartElement(reinterpret_cast<const char *>(localname), sPrefix, attributesMap);
}

void OnEndElementNs(
    void *         ctx,
    const xmlChar *localname,
    const xmlChar *prefix,
    const xmlChar *URI)
{
  ConfigParser *pParser = static_cast<ConfigParser *>(ctx);
  pParser->OnEndElement();
}

void OnCharacters(void *ctx, const xmlChar *ch, int len)
{
  ConfigParser *pParser = static_cast<ConfigParser *>(ctx);
  pParser->OnTextSection(std::string(reinterpret_cast<const char *>(ch), len));
}

void OnStructuredErrorFunc(void *userData, xmlError *error)
{
  const std::string message{error->message};

  // Ignore all namespace-related messages
  if (message.find("Namespace") != std::string::npos) {
    return;
  }

  ConfigParser::MessageProxy(error->level, message);
}

// ------------------------- ConfigParser implementation  -------------------------

precice::logging::Logger ConfigParser::_log("xml::XMLParser");

ConfigParser::ConfigParser(const std::string &filePath, const ConfigurationContext &context, std::shared_ptr<precice::xml::XMLTag> pXmlTag)
    : m_pXmlTag(std::move(pXmlTag))
{
  readXmlFile(filePath);

  std::vector<std::shared_ptr<XMLTag>> DefTags{m_pXmlTag};
  CTagPtrVec                           SubTags;
  // Initialize with the root tag, if any.
  if (not m_AllTags.empty())
    SubTags.push_back(m_AllTags[0]);

  try {
    connectTags(context, DefTags, SubTags);
  } catch (const std::exception &e) {
    PRECICE_ERROR("An unexpected exception occurred during configuration: " << e.what() << '.');
  }
}

ConfigParser::ConfigParser(const std::string &filePath)
{
  readXmlFile(filePath);
}

void ConfigParser::MessageProxy(int level, const std::string &mess)
{
  switch (level) {
  case (XML_ERR_FATAL):
  case (XML_ERR_ERROR):
    PRECICE_ERROR(mess);
    break;
  case (XML_ERR_WARNING):
    PRECICE_WARN(mess);
    break;
  default:
    PRECICE_INFO(mess);
  }
}

int ConfigParser::readXmlFile(std::string const &filePath)
{
  xmlSAXHandler SAXHandler;

  memset(&SAXHandler, 0, sizeof(xmlSAXHandler));

  SAXHandler.initialized    = XML_SAX2_MAGIC;
  SAXHandler.startElementNs = OnStartElementNs;
  SAXHandler.endElementNs   = OnEndElementNs;
  SAXHandler.characters     = OnCharacters;
  SAXHandler.serror         = OnStructuredErrorFunc;

  std::ifstream ifs(filePath);
  PRECICE_CHECK(ifs, "XML parser was unable to open configuration file \"" << filePath << '"');

  std::string content{std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>()};

  xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(&SAXHandler, static_cast<void *>(this),
                                                  content.c_str(), content.size(), nullptr);

  xmlParseChunk(ctxt, nullptr, 0, 1);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return 0;
}

void ConfigParser::connectTags(const ConfigurationContext &context, std::vector<std::shared_ptr<XMLTag>> &DefTags, CTagPtrVec &SubTags)
{
  std::unordered_set<std::string> usedTags;

  for (auto &subtag : SubTags) {
    std::string expectedName = (subtag->m_Prefix.length() ? subtag->m_Prefix + ":" : "") + subtag->m_Name;
    const auto  tagPosition  = std::find_if(
        DefTags.begin(),
        DefTags.end(),
        [expectedName](const std::shared_ptr<XMLTag> &pTag) {
          return pTag->_fullName == expectedName;
        });

    if (tagPosition == DefTags.end()) {
      PRECICE_ERROR("The configuration contains an unknown tag <" + expectedName + ">.");
    }

    auto pDefSubTag = *tagPosition;
    pDefSubTag->resetAttributes();

    if ((pDefSubTag->_occurrence == XMLTag::OCCUR_ONCE) || (pDefSubTag->_occurrence == XMLTag::OCCUR_NOT_OR_ONCE)) {
      if (usedTags.count(pDefSubTag->_fullName)) {
        PRECICE_ERROR("Tag <" + pDefSubTag->_fullName + "> is not allowed to occur multiple times.");
      }
      usedTags.emplace(pDefSubTag->_fullName);
    }

    pDefSubTag->_configuredNamespaces[pDefSubTag->_namespace] = true;
    pDefSubTag->readAttributes(subtag->m_aAttributes);
    pDefSubTag->_listener.xmlTagCallback(context, *pDefSubTag);
    pDefSubTag->_configured = true;

    connectTags(context, pDefSubTag->_subtags, subtag->m_aSubTags);

    pDefSubTag->areAllSubtagsConfigured();
    pDefSubTag->_listener.xmlEndTagCallback(context, *pDefSubTag);
  }
}

void ConfigParser::OnStartElement(
    std::string         localname,
    std::string         prefix,
    CTag::AttributePair attributes)
{
  auto pTag = std::make_shared<CTag>();

  pTag->m_Prefix      = std::move(prefix);
  pTag->m_Name        = std::move(localname);
  pTag->m_aAttributes = std::move(attributes);

  if (not m_CurrentTags.empty()) {
    auto pParentTag = m_CurrentTags.back();
    pParentTag->m_aSubTags.push_back(pTag);
  }

  m_AllTags.push_back(pTag);
  m_CurrentTags.push_back(pTag);
}

void ConfigParser::OnEndElement()
{
  m_CurrentTags.pop_back();
}

void ConfigParser::OnTextSection(const std::string &)
{
  // This page intentionally left blank
}
} // namespace xml
} // namespace precice

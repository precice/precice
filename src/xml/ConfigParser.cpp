#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iterator>
#include <libxml/SAX.h>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "utils/String.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLTag.hpp"

namespace precice::xml {

std::string decodeXML(std::string_view xml)
{
  static const std::map<std::string_view, char> escapes{{"&lt;", '<'}, {"&gt;", '>'}, {"&amp;", '&'}, {"&quot;", '"'}, {"&apos;", '\''}};
  std::string                                   decodedXml(xml);
  while (true) {
    bool changes{false};
    for (const auto &kv : escapes) {
      auto position = decodedXml.find(kv.first);
      if (position != std::string::npos) {
        decodedXml.replace(position, kv.first.length(), 1, kv.second);
        changes = true;
      }
    }
    if (!changes) {
      break;
    }
  };
  return decodedXml;
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

    auto             valueBegin = reinterpret_cast<const char *>(attributes[index + 3]);
    auto             valueEnd   = reinterpret_cast<const char *>(attributes[index + 4]);
    std::string_view value(valueBegin,
                           valueEnd - valueBegin);

    attributesMap[attributeName] = decodeXML(value);
  }

  auto pParser = static_cast<ConfigParser *>(ctx);

  std::string_view sPrefix(prefix == nullptr ? "" : reinterpret_cast<const char *>(prefix));

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

void OnStructuredErrorFunc(void *userData, const xmlError *error)
{
  const std::string_view message{error->message};

  // Ignore all namespace-related messages
  if (message.find("Namespace") != std::string::npos) {
    return;
  }

  ConfigParser::MessageProxy(error->level, message);
}

// Required for versions before 2.12.0 of libxml
void OnStructuredErrorFunc(void *userData, xmlError *error)
{
  OnStructuredErrorFunc(userData, static_cast<const xmlError *>(error));
}

void OnErrorFunc(void *userData, const char *error, ...)
{
  ConfigParser::MessageProxy(XML_ERR_ERROR, error);
}

void OnFatalErrorFunc(void *userData, const char *error, ...)
{
  ConfigParser::MessageProxy(XML_ERR_FATAL, error);
}

// ------------------------- ConfigParser implementation  -------------------------

precice::logging::Logger ConfigParser::_log("xml::XMLParser");

ConfigParser::ConfigParser(std::string_view filePath, const ConfigurationContext &context, std::shared_ptr<precice::xml::XMLTag> pXmlTag)
    : m_pXmlTag(std::move(pXmlTag))
{
  readXmlFile(std::string(filePath));

  std::vector<std::shared_ptr<XMLTag>> DefTags{m_pXmlTag};
  CTagPtrVec                           SubTags;
  // Initialize with the root tag, if any.
  if (not m_AllTags.empty())
    SubTags.push_back(m_AllTags[0]);

  try {
    connectTags(context, DefTags, SubTags);
  } catch (const std::exception &e) {
    PRECICE_ERROR("An unexpected exception occurred during configuration: {}.", e.what());
  }
}

ConfigParser::ConfigParser(std::string_view filePath)
{
  readXmlFile(std::string(filePath));
}

void ConfigParser::MessageProxy(int level, std::string_view mess)
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
  SAXHandler.error          = OnErrorFunc;
  SAXHandler.fatalError     = OnFatalErrorFunc;

  std::ifstream ifs{filePath};
  PRECICE_CHECK(ifs, "XML parser was unable to open configuration file \"{}\"", filePath);

  std::string content{std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>()};

  xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(&SAXHandler, static_cast<void *>(this),
                                                  content.c_str(), content.size(), nullptr);

  xmlParseChunk(ctxt, nullptr, 0, 1);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return 0;
}

namespace {
struct Distance {
  std::size_t distance;
  std::string name;

  bool operator<(const Distance &other) const
  {
    return distance < other.distance;
  }
};
auto gatherCandidates(const std::vector<std::shared_ptr<XMLTag>> &DefTags, std::string_view prefix)
{
  bool validPrefix = std::any_of(DefTags.begin(), DefTags.end(), [prefix](const auto &tag) { return tag->getNamespace() == prefix; });

  std::set<std::string> entries;
  for (const auto &tag : DefTags) {
    if (!validPrefix || (tag->getNamespace() == prefix)) {
      entries.insert(tag->getFullName());
    }
  }
  return entries;
}
} // namespace

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
      // Tag not found
      auto names = gatherCandidates(DefTags, subtag->m_Prefix);

      auto matches = utils::computeMatches(expectedName, names);
      if (!matches.empty() && matches.front().distance < 3) {
        matches.erase(std::remove_if(matches.begin(), matches.end(), [](auto &m) { return m.distance > 2; }), matches.end());
        std::vector<std::string_view> stringMatches;
        std::transform(matches.begin(), matches.end(), std::back_inserter(stringMatches), [](auto &m) { return m.name; });
        PRECICE_ERROR("The configuration contains an unknown tag <{}>. Did you mean <{}>?", expectedName, fmt::join(stringMatches, ">,<"));
      } else {
        PRECICE_ERROR("The configuration contains an unknown tag <{}>. Expected tags are {}.", expectedName, fmt::join(names, ", "));
      }
    }

    auto pDefSubTag = *tagPosition;
    pDefSubTag->resetAttributes();

    if ((pDefSubTag->_occurrence == XMLTag::OCCUR_ONCE) || (pDefSubTag->_occurrence == XMLTag::OCCUR_NOT_OR_ONCE)) {
      PRECICE_CHECK(usedTags.count(pDefSubTag->_fullName) == 0,
                    "Tag <{}> is not allowed to occur multiple times.", pDefSubTag->_fullName);
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
    std::string_view    localname,
    std::string_view    prefix,
    CTag::AttributePair attributes)
{
  auto pTag = std::make_shared<CTag>();

  pTag->m_Prefix      = prefix;
  pTag->m_Name        = localname;
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
} // namespace precice::xml

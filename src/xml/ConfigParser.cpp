#include "ConfigParser.hpp"
#include <libxml/SAX.h>
#include <fstream>

namespace precice
{
namespace xml
{

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

    const xmlChar *valueBegin = attributes[index + 3];
    const xmlChar *valueEnd   = attributes[index + 4];
    std::string    value((const char *) valueBegin, (const char *) valueEnd);

    attributesMap[attributeName] = value;
  }

  ConfigParser *pParser = static_cast<ConfigParser *>(ctx);

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

// ------------------------- ConfigParser implementation  -------------------------

precice::logging::Logger ConfigParser::_log("xml::XMLParser");

ConfigParser::ConfigParser(const std::string &filePath, std::shared_ptr<precice::xml::XMLTag> pXmlTag)
{
  m_pXmlTag = pXmlTag;
  readXmlFile(filePath);

  std::vector<std::shared_ptr<XMLTag>> DefTags;
  DefTags.push_back(m_pXmlTag);

  CTagPtrVec SubTags;
  // Initialize with the root tag, if any.
  if (not m_AllTags.empty())
    SubTags.push_back(m_AllTags[0]);

  try {
    connectTags(DefTags, SubTags);
  } catch (const std::string& error) {
    ERROR(error);
  }
}

ConfigParser::ConfigParser(const std::string &filePath)
{
  readXmlFile(filePath);
}

void ConfigParser::GenericErrorFunc(void *ctx, const char *msg, ...)
{
  const int TMP_BUF_SIZE = 256;

  char    err[TMP_BUF_SIZE];
  va_list arg_ptr;
  va_start(arg_ptr, msg);
  vsnprintf(err, TMP_BUF_SIZE, msg, arg_ptr);
  va_end(arg_ptr);

  ERROR(err);
}

int ConfigParser::readXmlFile(std::string const &filePath)
{
  auto handler = static_cast<xmlGenericErrorFunc>(ConfigParser::GenericErrorFunc);
  initGenericErrorDefaultFunc(&handler);

  xmlSAXHandler SAXHandler;

  memset(&SAXHandler, 0, sizeof(xmlSAXHandler));

  SAXHandler.initialized    = XML_SAX2_MAGIC;
  SAXHandler.startElementNs = OnStartElementNs;
  SAXHandler.endElementNs   = OnEndElementNs;
  SAXHandler.characters     = OnCharacters;

  std::ifstream ifs(filePath);
  if (not ifs) {
    ERROR("File open error: " << filePath);
  }

  std::string content{std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>()};

  xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(&SAXHandler, static_cast<void*>(this),
                                                  content.c_str(), content.size(), nullptr);

  xmlParseChunk(ctxt, nullptr, 0, 1);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return 0;
}

void ConfigParser::connectTags(std::vector<std::shared_ptr<XMLTag>> &DefTags, CTagPtrVec &SubTags)
{
  std::vector<std::string> usedTags;

  for (auto subtag : SubTags) {

    bool found = false;

    for (auto pDefSubTag : DefTags) {
      
      if (pDefSubTag->_fullName == ((subtag->m_Prefix.length() ? subtag->m_Prefix + ":" : "") + subtag->m_Name)) {
        found = true;
        pDefSubTag->resetAttributes();

        if (pDefSubTag->_occurrence == XMLTag::OCCUR_ONCE) {
          if (std::find(usedTags.begin(), usedTags.end(), pDefSubTag->_fullName) != usedTags.end()) {
            ERROR("Tag <" + pDefSubTag->_fullName + "> is already used");
          }
          usedTags.push_back(pDefSubTag->_fullName);
        }

        pDefSubTag->_configuredNamespaces[pDefSubTag->_namespace] = true;
        pDefSubTag->readAttributes(subtag->m_aAttributes);
        pDefSubTag->_listener.xmlTagCallback(*pDefSubTag);
        pDefSubTag->_configured = true;

        connectTags(pDefSubTag->_subtags, subtag->m_aSubTags);

        if (!pDefSubTag->_subtags.empty()) {
          pDefSubTag->areAllSubtagsConfigured();
          pDefSubTag->_listener.xmlEndTagCallback(*pDefSubTag);
        }

        break;
      }
    }

    if (!found)
      ERROR("Tag <" + ((subtag->m_Prefix.length() ? subtag->m_Prefix + ":" : "") + subtag->m_Name) + "> is unknown");
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

void ConfigParser::OnTextSection(std::string ch)
{
  // This page intentionally left blank
}
}
}

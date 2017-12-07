#include "ConfigParser.hpp"

namespace precice
{
namespace xml
{

precice::logging::Logger ConfigParser::_log("xml::XMLParser");

ConfigParser::ConfigParser(const std::string &filePath, precice::utils::XMLTag *pXmlTag)
{
  m_pXmlTag = pXmlTag;
  init(filePath);

  std::vector<precice::utils::XMLTag *> DefTags;
  DefTags.push_back(m_pXmlTag);

  std::vector<CTag *> SubTags;
  SubTags.push_back(getRootTag());

  try {
    connectTags(DefTags, SubTags);
  } catch (std::string error) {
    ERROR(error);
  }
}

ConfigParser::ConfigParser(const std::string &filePath)
{
  init(filePath);
}

void ConfigParser::init(const std::string &filePath)
{
  FILE *f = fopen(filePath.c_str(), "r");
  if (!f) {
    ERROR("file open error: " << filePath);
  }

  if (readXmlFile(f)) {
    ERROR("xml read error: " << filePath);
  }

  fclose(f);
}

ConfigParser::~ConfigParser()
{
  while (!m_AllTags.empty()) {
    CTag *pTag = m_AllTags.back();

    delete pTag;

    m_AllTags.pop_back();
  }
}

/*void GenericErrorFunc(void * ctx, const char * msg, ...)
{
	va_list args;
    va_start(args, msg);
 
    while (*msg != '\0') {
		std::string err = va_arg(args, const char *);
		std::cout << err;
		++msg;
	}
	
	va_end(args);
   
	std::cout << std::endl;
}*/

void ConfigParser::GenericErrorFunc(void *ctx, const char *msg, ...)
{
  const int TMP_BUF_SIZE = 256;

  char    err[TMP_BUF_SIZE];
  va_list arg_ptr;
  va_start(arg_ptr, msg);
  vsnprintf(err, TMP_BUF_SIZE, msg, arg_ptr);
  va_end(arg_ptr);

  ERROR("hm");
}

int ConfigParser::readXmlFile(FILE *f)
{
  char chars[1024];
  int  res = fread(chars, 1, 4, f);
  if (res <= 0) {
    return 1;
  }

  xmlGenericErrorFunc handler = (xmlGenericErrorFunc) ConfigParser::GenericErrorFunc;
  initGenericErrorDefaultFunc(&handler);

  xmlSAXHandler SAXHander = makeSaxHandler();

  xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(
      &SAXHander, (void *) this, chars, res, NULL);

  while ((res = fread(chars, 1, sizeof(chars), f)) > 0) {
    if (xmlParseChunk(ctxt, chars, res, 0)) {
      xmlParserError(ctxt, "xmlParseChunk");
      return 1;
    }
  }
  xmlParseChunk(ctxt, chars, 0, 1);

  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return 0;
}

ConfigParser::CTag *ConfigParser::getRootTag()
{
  if (m_AllTags.empty())
    return nullptr;

  return m_AllTags[0];
}

void ConfigParser::connectTags(std::vector<precice::utils::XMLTag *> &DefTags, std::vector<CTag *> &SubTags)
{
  std::vector<std::string> usedTags;

  for (auto subtag : SubTags) {
    CTag *pSubTag = (CTag *) subtag;

    bool found = false;

    for (auto defSubtag : DefTags) {
      precice::utils::XMLTag *pDefSubTag = defSubtag;

      if (pDefSubTag->_fullName == ((pSubTag->m_Prefix.length() ? pSubTag->m_Prefix + ":" : "") + pSubTag->m_Name)) {
        found = true;
        pDefSubTag->resetAttributes();

        if (pDefSubTag->_occurrence == precice::utils::XMLTag::OCCUR_ONCE) {
          if (std::find(usedTags.begin(), usedTags.end(), pDefSubTag->_fullName) != usedTags.end()) {
            ERROR("Tag <" + pDefSubTag->_fullName + "> is already used");
          }

          usedTags.push_back(pDefSubTag->_fullName);
        }

        pDefSubTag->_configuredNamespaces[pDefSubTag->_namespace] = true;
        pDefSubTag->readAttributes(pSubTag->m_aAttributes);
        pDefSubTag->_listener.xmlTagCallback(*pDefSubTag);
        pDefSubTag->_configured = true;

        connectTags(pDefSubTag->_subtags, pSubTag->m_aSubTags);

        if (!pDefSubTag->_subtags.empty()) {
          pDefSubTag->areAllSubtagsConfigured();
          pDefSubTag->_listener.xmlEndTagCallback(*pDefSubTag);
        }

        break;
      }
    }

    if (!found)
      ERROR("Tag <" + ((pSubTag->m_Prefix.length() ? pSubTag->m_Prefix + ":" : "") + pSubTag->m_Name) + "> is unknown");
  }
}

xmlSAXHandler ConfigParser::makeSaxHandler()
{
  xmlSAXHandler SAXHander;

  memset(&SAXHander, 0, sizeof(xmlSAXHandler));

  SAXHander.initialized    = XML_SAX2_MAGIC;
  SAXHander.startElementNs = OnStartElementNs;
  SAXHander.endElementNs   = OnEndElementNs;
  SAXHander.characters     = OnCharacters;

  return SAXHander;
}

void ConfigParser::OnStartElementNs(
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
  CTag *pTag = new CTag();

  if (prefix != nullptr)
    pTag->m_Prefix = std::string((const char *) prefix);

  pTag->m_Name = std::string((const char *) localname);

  unsigned int index = 0;
  for (int indexAttribute = 0; indexAttribute < nb_attributes; ++indexAttribute, index += 5) {
    const xmlChar *localname  = attributes[index];
    const xmlChar *valueBegin = attributes[index + 3];
    const xmlChar *valueEnd   = attributes[index + 4];

    std::string value((const char *) valueBegin, (const char *) valueEnd);

    pTag->m_aAttributes[std::string((const char *) localname)] = value;
  }

  ConfigParser *pParser = (ConfigParser *) ctx;

  if (!pParser->m_CurrentTags.empty()) {
    CTag *pParentTag = pParser->m_CurrentTags.back();
    pParentTag->m_aSubTags.push_back(pTag);
  }

  pParser->m_AllTags.push_back(pTag);
  pParser->m_CurrentTags.push_back(pTag);
}

void ConfigParser::OnEndElementNs(
    void *         ctx,
    const xmlChar *localname,
    const xmlChar *prefix,
    const xmlChar *URI)
{
  ConfigParser *pParser = (ConfigParser *) ctx;

  pParser->m_CurrentTags.pop_back();
}

void ConfigParser::OnCharacters(void *ctx, const xmlChar *ch, int len)
{
  char chars[len + 1];
  strncpy(chars, (const char *) ch, len);
  chars[len] = (char) NULL;
  //std::cout << "Text: " << chars << std::endl;
}
}
}

#pragma once

#include <libxml/parser.h>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace logging {
class Logger;
} // namespace logging

namespace xml {
class XMLTag; // forward declaration to resolve circular import
struct ConfigurationContext;

/// Decodes escape sequences of a given xml
std::string decodeXML(std::string_view xml);

class ConfigParser {
public:
  /// Struct holding the read tag from xml file
  struct CTag {
    std::string m_Name;
    std::string m_Prefix;
    bool        m_Used = false;

    using AttributePair = std::map<std::string, std::string>;
    AttributePair                      m_aAttributes;
    std::vector<std::shared_ptr<CTag>> m_aSubTags;

    /// 1-based line and column of the opening tag in the configuration file, or -1 if unknown.
    int m_Line   = -1;
    int m_Column = -1;
  };

  using CTagPtrVec = std::vector<std::shared_ptr<CTag>>;

private:
  static precice::logging::Logger _log;

  /// the hash of the last processed config
  std::string _hash;

  CTagPtrVec m_AllTags;
  CTagPtrVec m_CurrentTags;

  std::shared_ptr<precice::xml::XMLTag> m_pXmlTag;

  /// Lines of the configuration file currently being parsed, used to add source context to error messages.
  std::vector<std::string> m_lines;

  /// The libxml2 parser context, only valid while readXmlFile() is parsing.
  xmlParserCtxtPtr m_parserCtxt = nullptr;

  /// Returns the source line (1-based), or an empty string if unknown/out-of-range.
  std::string_view sourceLine(int line) const;

public:
  /// Parser ctor for Callback init
  ConfigParser(std::string_view filePath, const ConfigurationContext &context, std::shared_ptr<XMLTag> pXmlTag);

  /// Parser ctor without Callbacks
  ConfigParser(std::string_view filePath);

  /// Reads the xml file
  int readXmlFile(std::string const &filePath);

  /// Reads and returns the content of the file.
  std::string readFileContent(std::string const &filePath) const;

  /// returns the hash of the processed XML file
  std::string hash() const;

  /**
   * @brief Connects the actual tags of an xml layer with the predefined tags
   * @param DefTags predefined tags
   * @param SubTags actual tags from xml file
   */
  void connectTags(const ConfigurationContext &context, std::vector<std::shared_ptr<precice::xml::XMLTag>> &DefTags, CTagPtrVec &SubTags);

  /// Callback for Start-Tag
  void OnStartElement(
      std::string_view    localname,
      std::string_view    prefix,
      CTag::AttributePair attributes,
      int                 line,
      int                 column);

  /// Callback for End-Tag
  void OnEndElement();

  /// Callback for text sections in xml file
  void OnTextSection(const std::string &ch);

  /// Returns the current line/column of the parser, or {-1, -1} if not currently parsing.
  std::pair<int, int> currentParserPosition() const;

  /// Reports an error or warning message from libxml2, annotated with the given source location.
  void reportParserMessage(int level, std::string_view message, int line, int column) const;
};
} // namespace xml
} // namespace precice

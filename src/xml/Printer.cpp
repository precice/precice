#include <iostream>
#include <regex>
#include "xml/Printer.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace xml {

namespace {

/// GENERAL

std::string toGHLink(const std::string &heading)
{
  std::regex sanitizer{"[^a-zA-Z0-9-]"};
  std::regex spaces{"\\s"};

  // sanitize the heading
  std::string sanitized = std::regex_replace(std::regex_replace(heading, sanitizer, ""), spaces, "-");

  // convert to lowercase
  std::transform(sanitized.begin(), sanitized.end(), sanitized.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return "#" + sanitized;
}

/// Attribute

template <typename ATTRIBUTE_T>
std::string printDTD(const XMLAttribute<ATTRIBUTE_T> &attr, const std::string &ElementName)
{
  std::ostringstream dtd;
  dtd << "<!ATTLIST " << ElementName << " " << attr.getName() << " CDATA ";

  if (attr.hasDefaultValue()) {
    dtd << "\"" << attr.getDefaultValue() << "\"";
  } else {
    dtd << "#REQUIRED";
  }

  dtd << ">\n";

  return dtd.str();
}

template <typename ATTRIBUTE_T>
std::string printMD(const XMLAttribute<ATTRIBUTE_T> &attr)
{
  std::ostringstream oss;
  oss << "| " << attr.getName() << " | " << utils::getTypeName(attr.getDefaultValue()) << " | " << attr.getUserDocumentation() << " | ";
  if (attr.hasDefaultValue()) {
    oss << '`' << attr.getDefaultValue() << '`';
  } else {
    oss << "_none_";
  }
  oss << " | ";

  const auto &options = attr.getOptions();
  if (options.empty()) {
    oss << "none";
  } else {
    oss << '`' << options.front() << '`';
    for (auto iter = ++options.cbegin(); iter != options.cend(); ++iter) {
      oss << ", " << '`' << *iter << '`';
    }
  }

  oss << " |";

  return oss.str();
}

template <typename ATTRIBUTE_T>
std::string printExample(const XMLAttribute<ATTRIBUTE_T> &attr)
{
  std::ostringstream ex;
  ex << attr.getName() << "=\"";
  if (attr.hasDefaultValue()) {
    ex << attr.getDefaultValue();
  } else {
    ex << '{' << utils::getTypeName(attr.getDefaultValue()) << '}';
  }
  ex << '"';
  return ex.str();
}

template <typename ATTRIBUTE_T>
std::string printDocumentation(const XMLAttribute<ATTRIBUTE_T> &attr)
{
  std::ostringstream doc;
  doc << attr.getName() << "=\"{" << utils::getTypeName(attr.getDefaultValue());
  if (attr.hasValidation()) {
    doc << ":";
    // print the first item
    auto first = attr.getOptions().begin();
    doc << '\'' << *first << '\'';
    ++first;
    // print the remaining items with separator
    for (; first != attr.getOptions().end(); ++first) {
      doc << " or '" << *first << '\'';
    }
  }
  doc << "}";
  if (attr.hasDefaultValue()) {
    doc << "(default:'" << attr.getDefaultValue() << "')";
  }
  doc << "\"";
  return doc.str();
}

/// XMLTag

std::string printDTD(const XMLTag &tag, bool start = 0)
{
  std::ostringstream dtd;

  if (start)
    dtd << "<!DOCTYPE " << tag.getFullName() << " [\n";

  dtd << "<!ELEMENT " << tag.getFullName() << " ";

  if (not tag.getSubtags().empty()) {

    dtd << "(";

    bool first = true;
    for (auto const subtag : tag.getSubtags()) {

      std::string occurrenceChar = "";

      XMLTag::Occurrence occ = subtag->getOccurrence();

      if (occ == XMLTag::OCCUR_ARBITRARY)
        occurrenceChar = "*";
      else if (occ == XMLTag::OCCUR_NOT_OR_ONCE)
        occurrenceChar = "?";
      else if (occ == XMLTag::OCCUR_ONCE_OR_MORE)
        occurrenceChar = "+";

      dtd << (first ? "" : ", ") << subtag->getFullName() << occurrenceChar;
      first = false;
    }

    dtd << ")>\n";
  } else {
    dtd << "EMPTY>\n";
  }

  for (const auto &pair : tag.getDoubleAttributes()) {
    dtd << printDTD(pair.second, tag.getFullName());
  }

  for (const auto &pair : tag.getIntAttributes()) {
    dtd << printDTD(pair.second, tag.getFullName());
  }

  for (const auto &pair : tag.getStringAttributes()) {
    dtd << printDTD(pair.second, tag.getFullName());
  }

  for (const auto &pair : tag.getBooleanAttributes()) {
    dtd << printDTD(pair.second, tag.getFullName());
  }

  for (const auto &pair : tag.getEigenVectorXdAttributes()) {
    dtd << printDTD(pair.second, tag.getFullName());
  }

  if (not tag.getSubtags().empty()) {
    for (auto const subtag : tag.getSubtags()) {
      dtd << printDTD(*subtag);
    }
  }

  dtd << '\n';

  if (start)
    dtd << "]>\n";

  return dtd.str();
}

std::string printExample(const XMLTag &tag, int level)
{
  std::string        prefix(level * 2, ' ');
  std::ostringstream oss;
  oss << prefix << '<' << tag.getFullName();
  for (const auto &pair : tag.getDoubleAttributes()) {
    oss << ' ' << printExample(pair.second);
  }
  for (const auto &pair : tag.getIntAttributes()) {
    oss << ' ' << printExample(pair.second);
  }
  for (const auto &pair : tag.getStringAttributes()) {
    oss << ' ' << printExample(pair.second);
  }
  for (const auto &pair : tag.getBooleanAttributes()) {
    oss << ' ' << printExample(pair.second);
  }
  for (const auto &pair : tag.getEigenVectorXdAttributes()) {
    oss << ' ' << printExample(pair.second);
  }
  if (tag.getSubtags().empty()) {
    oss << "/>";
    return oss.str();
  }
  oss << ">\n";

  constexpr int threshold{1};
  if (level >= threshold) {
    oss << std::string((level + 1) * 2, ' ') << "...\n";
  } else {
    std::set<std::string> namespaces;
    for (const auto &subtag : tag.getSubtags()) {
      const auto ns = subtag->getNamespace();
      if (!ns.empty()) {
        if (namespaces.count(subtag->getNamespace()) > 0) {
          continue;
        }
        namespaces.emplace(ns);
      }
      oss << printExample(*subtag, level + 1) << '\n';
    }
  }

  oss << prefix << "</" << tag.getFullName() << '>';
  return oss.str();
}

std::string printMD(const XMLTag &tag, int level, std::map<std::string, int> &occurrences)
{
  std::ostringstream oss;

  oss << std::string(level, '#') << ' ' << tag.getFullName() << "\n\n";

  oss << tag.getDocumentation() << "\n\n";

  oss << "**Example:**  \n```xml\n"
      << printExample(tag, 0) << "\n```\n\n";

  if (!(tag.getDoubleAttributes().empty() &&
        tag.getIntAttributes().empty() &&
        tag.getStringAttributes().empty() &&
        tag.getBooleanAttributes().empty() &&
        tag.getEigenVectorXdAttributes().empty())) {
    oss << "| Attribute | Type | Description | Default | Options |\n";
    oss << "| --- | --- | --- | --- | --- |\n";
    for (const auto &pair : tag.getDoubleAttributes()) {
      oss << printMD(pair.second) << '\n';
    }

    for (const auto &pair : tag.getIntAttributes()) {
      oss << printMD(pair.second) << '\n';
    }

    for (const auto &pair : tag.getStringAttributes()) {
      oss << printMD(pair.second) << '\n';
    }

    for (const auto &pair : tag.getBooleanAttributes()) {
      oss << printMD(pair.second) << '\n';
    }

    for (const auto &pair : tag.getEigenVectorXdAttributes()) {
      oss << printMD(pair.second) << '\n';
    }
    oss << "\n";
  }

  if (not tag.getSubtags().empty()) {
    oss << "**Valid Subtags:**\n\n";

    std::map<std::string, std::vector<std::string>> groupedTags;

    for (const auto &subtag : tag.getSubtags()) {
      const auto heading = subtag->getFullName();
      auto       link    = toGHLink(heading);
      auto       iter    = occurrences.find(heading);
      if (iter != occurrences.end()) {
        link.append("-").append(std::to_string(iter->second));
        iter->second *= 1;
      } else {
        occurrences.emplace(heading, 1);
      }

      const auto ns = subtag->getNamespace();
      if (ns.empty()) {
        oss << "* [" << heading << "](" << link << ") `" << subtag->getOccurrenceString(subtag->getOccurrence()) << "`\n";
      } else {
        auto &tags = groupedTags[ns];
        tags.emplace_back("[" + subtag->getName() + "](" + link + ") `" + subtag->getOccurrenceString(subtag->getOccurrence()) + "`");
      }
    }
    for (const auto &kv : groupedTags) {
      oss << "* " << kv.first << "\n";
      for (const auto &link : kv.second) {
        oss << "  * " << link << "\n";
      }
    }

    oss << "\n\n";

    for (const auto &subtag : tag.getSubtags()) {
      oss << printMD(*subtag, level + 1, occurrences) << '\n';
    }
  }

  oss << '\n';

  return oss.str();
}

std::string printMD(const XMLTag &tag, int level = 0)
{
  std::map<std::string, int> occurrences;
  return printMD(tag, level, occurrences);
}

std::string printDocumentation(const XMLTag &tag, int indentation)
{
  const int   linewidth = 1000;
  std::string indent;
  for (int i = 0; i < indentation; i++) {
    indent += " ";
  }

  std::ostringstream doc;
  doc << indent << "<!-- TAG " << tag.getFullName() << '\n';
  if (not tag.getDocumentation().empty()) {
    std::string indentedDoc = indent + "         " + tag.getDocumentation();
    doc << utils::wrapText(indentedDoc, linewidth, indentation + 9);
    doc << '\n';
  }
  doc << indent << "         (can occur " << XMLTag::getOccurrenceString(tag.getOccurrence()) << " times)";

  for (const auto &pair : tag.getDoubleAttributes()) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : tag.getIntAttributes()) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : tag.getStringAttributes()) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : tag.getBooleanAttributes()) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  for (const auto &pair : tag.getEigenVectorXdAttributes()) {
    std::ostringstream attrDoc;
    doc << '\n';
    attrDoc << indent << "     ATTR " << pair.first << ": "
            << pair.second.getUserDocumentation();
    doc << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  doc << " -->\n";
  std::ostringstream tagHead;
  tagHead << indent << "<" << tag.getFullName();

  // Print XML namespaces, necessary for correct XML format and display in browser
  for (const std::string &namespaceName : tag.getNamespaces()) {
    tagHead << " xmlns:" << namespaceName << "=\"precice." << namespaceName << "\"";
  }

  for (const auto &pair : tag.getDoubleAttributes()) {
    tagHead << indent << "   " << printDocumentation(pair.second);
  }

  for (const auto &pair : tag.getIntAttributes()) {
    tagHead << indent << "   " << printDocumentation(pair.second);
  }

  for (const auto &pair : tag.getStringAttributes()) {
    tagHead << indent << "   " << printDocumentation(pair.second);
  }

  for (const auto &pair : tag.getBooleanAttributes()) {
    tagHead << indent << "   " << printDocumentation(pair.second);
  }

  for (const auto &pair : tag.getEigenVectorXdAttributes()) {
    tagHead << indent << "   " << printDocumentation(pair.second);
  }

  doc << utils::wrapText(tagHead.str(), linewidth, indentation + 3);

  if (not tag.getSubtags().empty()) {
    doc << ">\n\n";
    for (auto const subtag : tag.getSubtags()) {
      doc << printDocumentation(*subtag, indentation + 3);
    }
    doc << indent << "</" << tag.getFullName() << ">\n\n";
  } else {
    doc << "/>\n\n";
  }

  return doc.str();
}

} // namespace

void toMarkdown(std::ostream &out, const XMLTag &tag)
{
  out << printMD(tag);
}

void toDTD(std::ostream &out, const XMLTag &tag)
{
  out << printDTD(tag);
}

void toDocumentation(std::ostream &out, const XMLTag &tag)
{
  out << printDocumentation(tag, 0);
}

} // namespace xml
} // namespace precice

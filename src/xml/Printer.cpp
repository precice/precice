#include "xml/Printer.hpp"
#include <algorithm>
#include <cctype>
#include <fmt/format.h>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "utils/String.hpp"
#include "utils/TypeNames.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice::xml {

namespace {

/// GENERAL

/// Transforms a given heading into the id of the sanitized GitHub link
std::string toGHLink(const std::string &heading)
{
  try {
    std::regex sanitizer{"[^a-zA-Z0-9-]"};
    std::regex spaces{"\\s"};

    // sanitize the heading
    std::string sanitized = std::regex_replace(std::regex_replace(heading, sanitizer, ""), spaces, "-");

    // convert to lowercase
    std::transform(sanitized.begin(), sanitized.end(), sanitized.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return "#" + sanitized;

  } catch (const std::regex_error &e) {
    std::cerr << "Error sanitizing link: " << e.what() << '\n';
    std::exit(-1);
  }
}

/// Attribute

/// Prints the DTD attribute spec for an XMLAttribute given its name.
template <typename ATTRIBUTE_T>
std::ostream &printDTD(std::ostream &out, const XMLAttribute<ATTRIBUTE_T> &attr, const std::string &ElementName)
{
  out << "<!ATTLIST " << ElementName << " " << attr.getName() << " CDATA ";

  if (attr.hasDefaultValue()) {
    out << "\"" << attr.getDefaultValue() << "\"";
  } else {
    out << "#REQUIRED";
  }

  out << ">\n";
  return out;
}

/// Prints the attribute as a markdown table row.
template <typename ATTRIBUTE_T>
std::ostream &printMD(std::ostream &out, const XMLAttribute<ATTRIBUTE_T> &attr)
{
  fmt::print(out,
             "| {} | {} | {} |",
             attr.getName(),
             utils::getTypeName(attr.getDefaultValue()),
             attr.getUserDocumentation());

  if (attr.hasDefaultValue()) {
    fmt::print(out, " `{}` |", attr.getDefaultValue());
  } else {
    out << " _none_ |";
  }

  const auto &options = attr.getOptions();
  if (options.empty()) {
    out << " none |";
  } else {
    fmt::print(out, " `{}` |", fmt::join(options, "`, `"));
  }
  return out;
}

/// Prints an example of a given XMLAttribute.
template <typename ATTRIBUTE_T>
std::ostream &printExample(std::ostream &out, const XMLAttribute<ATTRIBUTE_T> &attr)
{
  out << attr.getName() << "=\"";
  if (attr.hasDefaultValue()) {
    fmt::print(out, "{}", attr.getDefaultValue());
  } else {
    out << '{' << utils::getTypeName(attr.getDefaultValue()) << '}';
  }
  out << '"';
  return out;
}

/// Prints the xml documentation of a given XMLAttribute
template <typename ATTRIBUTE_T>
std::ostream &printDocumentation(std::ostream &out, const XMLAttribute<ATTRIBUTE_T> &attr)
{
  out << attr.getName() << "=\"{" << utils::getTypeName(attr.getDefaultValue());
  if (attr.hasValidation()) {
    out << ":";
    // print the first item
    auto first = attr.getOptions().begin();
    out << '\'' << *first << '\'';
    ++first;
    // print the remaining items with separator
    for (; first != attr.getOptions().end(); ++first) {
      out << " or '" << *first << '\'';
    }
  }
  out << "}";
  if (attr.hasDefaultValue()) {
    out << "(default:'" << attr.getDefaultValue() << "')";
  }
  out << "\"";
  return out;
}

/// XMLTag

/** Prints the dtd of a given XMLTag.
 * Also prints the doctype if start is true.
 */
std::ostream &printDTD(std::ostream &out, const XMLTag &tag, bool start = false)
{
  if (start)
    out << "<!DOCTYPE " << tag.getFullName() << " [\n";

  out << "<!ELEMENT " << tag.getFullName() << " ";

  if (not tag.getSubtags().empty()) {

    out << "(";

    bool first = true;
    for (auto const &subtag : tag.getSubtags()) {

      std::string occurrenceChar;

      XMLTag::Occurrence occ = subtag->getOccurrence();

      if (occ == XMLTag::OCCUR_ARBITRARY)
        occurrenceChar = "*";
      else if (occ == XMLTag::OCCUR_NOT_OR_ONCE)
        occurrenceChar = "?";
      else if (occ == XMLTag::OCCUR_ONCE_OR_MORE)
        occurrenceChar = "+";

      out << (first ? "" : ", ") << subtag->getFullName() << occurrenceChar;
      first = false;
    }

    out << ")>\n";
  } else {
    out << "EMPTY>\n";
  }

  for (const auto &attribute : tag.getAttributes()) {
    std::visit([&out, &tag](const auto &attribute) { printDTD(out, attribute, tag.getFullName()); }, attribute);
  }

  if (not tag.getSubtags().empty()) {
    for (const auto &subtag : tag.getSubtags()) {
      printDTD(out, *subtag);
    }
  }

  out << '\n';

  if (start)
    out << "]>\n";
  return out;
}

/** Prints an Example of the given XMLTag at a given level of nesting.
 *
 * For the sake of readability, the example depth is truncated.
 * level is used to truncate and for indentation purposes.
 */
std::ostream &printExample(std::ostream &out, const XMLTag &tag, int level)
{
  std::string prefix(level * 2, ' ');
  out << prefix << '<' << tag.getFullName();
  for (const auto &attribute : tag.getAttributes()) {
    out << ' ';
    std::visit([&out](const auto &attribute) { printExample(out, attribute); }, attribute);
  }
  if (tag.getSubtags().empty()) {
    out << "/>";
    return out;
  }
  out << ">\n";

  constexpr int threshold{1};
  if (level >= threshold) {
    out << std::string((level + 1) * 2, ' ') << "...\n";
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
      printExample(out, *subtag, level + 1) << '\n';
    }
  }

  out << prefix << "</" << tag.getFullName() << '>';
  return out;
}

/** Prints the markdown reference for a given XMLTag
 *
 * @param[in] tag the tag to print
 * @param[in] level the level of nesting
 * @param[in] occurrences the occuences of tags, required to generate links
 */
std::ostream &printMD(std::ostream &out, const XMLTag &tag, int level, std::map<std::string, int> &occurrences)
{
  out << std::string(level, '#') << ' ' << tag.getFullName() << "\n\n";

  out << tag.getDocumentation() << "\n\n";

  out << "**Example:**  \n```xml\n";
  printExample(out, tag, 0) << "\n```\n\n";

  if (const auto &attributes = tag.getAttributes();
      !attributes.empty()) {
    out << "| Attribute | Type | Description | Default | Options |\n";
    out << "| --- | --- | --- | --- | --- |\n";

    for (const auto &attribute : attributes) {
      std::visit([&out](const auto &attribute) { printMD(out, attribute) << '\n'; }, attribute);
    }
    out << '\n';
  }

  if (not tag.getSubtags().empty()) {
    out << "**Valid Subtags:**\n\n";

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
        fmt::print(out, "* [{}]({}) `{}`\n",
                   heading,
                   link,
                   subtag->getOccurrenceString(subtag->getOccurrence()));
      } else {
        auto &tags = groupedTags[ns];
        tags.emplace_back(fmt::format("[{}]({}) `{}`",
                                      subtag->getName(),
                                      link,
                                      subtag->getOccurrenceString(subtag->getOccurrence())));
      }
    }
    for (const auto &kv : groupedTags) {
      out << "* " << kv.first << "\n";
      for (const auto &link : kv.second) {
        out << "  * " << link << "\n";
      }
    }

    out << "\n\n";

    for (const auto &subtag : tag.getSubtags()) {
      printMD(out, *subtag, level + 1, occurrences) << '\n';
    }
  }

  out << '\n';
  return out;
}

/** Print the markdown reference for a given XMLTag.
 * @note Use this as an initial call
 */
std::ostream &printMD(std::ostream &out, const XMLTag &tag, int level = 1)
{
  std::map<std::string, int> occurrences;
  printMD(out, tag, level, occurrences);
  return out;
}

/// Print the xml documentation of an XMLTag at a given level of indentation.
std::ostream &printDocumentation(std::ostream &out, const XMLTag &tag, int indentation)
{
  const int   linewidth = 1000;
  std::string indent;
  for (int i = 0; i < indentation; i++) {
    indent += " ";
  }

  out << indent << "<!-- TAG " << tag.getFullName() << '\n';
  if (not tag.getDocumentation().empty()) {
    std::string indentedDoc = indent + "         " + tag.getDocumentation();
    out << utils::wrapText(indentedDoc, linewidth, indentation + 9);
    out << '\n';
  }
  out << indent << "         (can occur " << XMLTag::getOccurrenceString(tag.getOccurrence()) << " times)";

  for (const auto &attribute : tag.getAttributes()) {
    out << '\n';
    std::ostringstream attrDoc;
    attrDoc << indent << "     ATTR " << getName(attribute) << ": "
            << std::visit([](const auto &attribute) { return attribute.getUserDocumentation(); }, attribute);
    out << utils::wrapText(attrDoc.str(), linewidth, indentation + 10);
  }

  out << " -->\n";
  std::ostringstream tagHead;
  tagHead << indent << "<" << tag.getFullName();

  // Print XML namespaces, necessary for correct XML format and display in browser
  for (const std::string &namespaceName : tag.getNamespaces()) {
    tagHead << " xmlns:" << namespaceName << "=\"precice." << namespaceName << "\"";
  }

  for (const auto &attribute : tag.getAttributes()) {
    tagHead << indent << "   ";
    std::visit([&tagHead](const auto &attribute) { printDocumentation(tagHead, attribute); }, attribute);
  }

  out << utils::wrapText(tagHead.str(), linewidth, indentation + 3);

  if (not tag.getSubtags().empty()) {
    out << ">\n\n";
    for (const auto &subtag : tag.getSubtags()) {
      printDocumentation(out, *subtag, indentation + 3);
    }
    out << indent << "</" << tag.getFullName() << ">\n\n";
  } else {
    out << "/>\n\n";
  }
  return out;
}

} // namespace

void toMarkdown(std::ostream &out, const XMLTag &tag)
{
  printMD(out, tag);
}

void toDTD(std::ostream &out, const XMLTag &tag)
{
  printDTD(out, tag);
}

void toDocumentation(std::ostream &out, const XMLTag &tag)
{
  printDocumentation(out, tag, 0);
}

} // namespace precice::xml

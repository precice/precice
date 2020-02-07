#pragma once

#include <iosfwd>

namespace precice {
namespace xml {

class XMLTag;

/// Prints the Markdown reference for the given tag.
void toMarkdown(std::ostream &out, const XMLTag &tag);

/// Prints the DTD reference for the given tag.
void toDTD(std::ostream &out, const XMLTag &tag);

/// Prints the XML reference for the given tag.
void toDocumentation(std::ostream &out, const XMLTag &tag);

} // namespace xml
} // namespace precice

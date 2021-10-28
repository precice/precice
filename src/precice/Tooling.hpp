#pragma once

#include <iosfwd>
#include <string>

namespace precice {

namespace tooling {

/// @brief Prints the configuration reference as Markdown.
void printConfigAsMD(std::ostream &out);

/// @brief Prints the configuration reference as DTD.
void printConfigAsDTD(std::ostream &out);

/// @brief Prints the configuration reference as XML with inlined help.
void printConfigAsXML(std::ostream &out);

} // namespace tooling

} // namespace precice

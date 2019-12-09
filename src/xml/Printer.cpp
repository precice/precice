#include <iostream>
#include "xml/Printer.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace xml {

void toMarkdown(std::ostream &out, const XMLTag &tag)
{
  out << tag.printMD();
}

void toDTD(std::ostream &out, const XMLTag &tag)
{
  out << tag.printDTD();
}

void toDocumentation(std::ostream &out, const XMLTag &tag)
{
  out << tag.printDocumentation(0);
}

} // namespace xml
} // namespace precice

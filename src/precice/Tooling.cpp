#include "precice/Tooling.hpp"

#include "precice/config/Configuration.hpp"
#include "precice/impl/versions.hpp"
#include "xml/Printer.hpp"

namespace precice {

namespace tooling {

void printConfigAsMD(std::ostream &out)
{
  precice::config::Configuration config;
  out << "<!-- generated with preCICE " PRECICE_VERSION " -->\n";
  precice::xml::toMarkdown(out, config.getXMLTag());
}

void printConfigAsDTD(std::ostream &out)
{
  precice::config::Configuration config;
  precice::xml::toDTD(out, config.getXMLTag());
}

void printConfigAsXML(std::ostream &out)
{
  precice::config::Configuration config;
  precice::xml::toDocumentation(out, config.getXMLTag());
}

} // namespace tooling

} // namespace precice

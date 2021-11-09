#include "precice/Tooling.hpp"

#include "precice/config/Configuration.hpp"
#include "precice/impl/versions.hpp"
#include "xml/Printer.hpp"

namespace precice {

namespace tooling {

void printConfigReference(std::ostream &out, ConfigReferenceType reftype)
{
  precice::config::Configuration config;
  switch (reftype) {
  case ConfigReferenceType::XML:
    precice::xml::toDocumentation(out, config.getXMLTag());
    return;
  case ConfigReferenceType::DTD:
    precice::xml::toDTD(out, config.getXMLTag());
    return;
  case ConfigReferenceType::MD:
    out << "<!-- generated with preCICE " PRECICE_VERSION " -->\n";
    precice::xml::toMarkdown(out, config.getXMLTag());
    return;
  }
}

} // namespace tooling

} // namespace precice

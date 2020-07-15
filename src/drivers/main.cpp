#include <iostream>
#include <string>
#include "precice/config/Configuration.hpp"
#include "utils/assertion.hpp"
#include "xml/Printer.hpp"

void printUsage()
{
  std::cout << "Usage:\n\n";
  std::cout << "Print XML reference      :  binprecice xml\n";
  std::cout << "Print DTD for XML config :  binprecice dtd" << std::endl;
  std::cout << "Print Markdown reference :  binprecice md" << std::endl;
}

int main(int argc, char **argv)
{
  bool runHelp = false;
  bool runDtd  = false;
  bool runMD   = false;

  bool wrongParameters = true;

  if (argc >= 2) {
    std::string action(argv[1]);
    if (action == "dtd") {
      wrongParameters = false;
      runDtd          = true;
    }
    if (action == "md") {
      wrongParameters = false;
      runMD           = true;
    }
    if (action == "xml") {
      wrongParameters = false;
      runHelp         = true;
    }
  }

  if (wrongParameters) {
    printUsage();
    return 1;
  }

  if (runHelp) {
    precice::config::Configuration config;
    precice::xml::toDocumentation(std::cout, config.getXMLTag());
  } else if (runDtd) {
    precice::config::Configuration config;
    precice::xml::toDTD(std::cout, config.getXMLTag());
  } else if (runMD) {
    precice::config::Configuration config;
    precice::xml::toMarkdown(std::cout, config.getXMLTag());
  } else {
    PRECICE_ASSERT(false);
  }
  return 0;
}

#include <cassert>
#include <iostream>
#include <precice/SolverInterface.hpp>
#include <string>

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
    precice::printConfigAsXML(std::cout);
  } else if (runDtd) {
    precice::printConfigAsDTD(std::cout);
  } else if (runMD) {
    precice::printConfigAsMD(std::cout);
  } else {
    assert(false);
  }
  return 0;
}

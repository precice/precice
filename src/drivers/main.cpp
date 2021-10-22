#include <cassert>
#include <iostream>
#include <precice/SolverInterface.hpp>
#include <string>

void printUsage()
{
  std::cerr << "Usage:\n\n";
  std::cerr << "Print XML reference      :  binprecice xml\n";
  std::cerr << "Print DTD for XML config :  binprecice dtd\n";
  std::cerr << "Print Markdown reference :  binprecice md\n";
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    printUsage();
    return 1;
  }

  const std::string action(argv[1]);
  const int         args = argc - 2;

  if (action == "dtd" && args == 0) {
    precice::printConfigAsXML(std::cout);
    return 0;
  }
  if (action == "md" && args == 0) {
    precice::printConfigAsMD(std::cout);
    return 0;
  }
  if (action == "xml" && args == 0) {
    precice::printConfigAsDTD(std::cout);
    return 0;
  }

  printUsage();
  return 1;
}

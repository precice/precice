#include <cassert>
#include <iostream>
#include <precice/SolverInterface.hpp>
#include <precice/Tooling.hpp>
#include <stdexcept>
#include <string>

void printUsage()
{
  std::cerr << "Usage:\n\n";
  std::cerr << "Print XML reference      :  binprecice xml\n";
  std::cerr << "Print DTD for XML config :  binprecice dtd\n";
  std::cerr << "Print Markdown reference :  binprecice md\n";
  std::cerr << "Print preCICE version    :  binprecice version\n";
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    printUsage();
    return 1;
  }

  const std::string action(argv[1]);
  const int         args = argc - 2;

  using namespace precice::tooling;

  if (action == "dtd" && args == 0) {
    printConfigReference(std::cout, ConfigReferenceType::DTD);
    return 0;
  }
  if (action == "md" && args == 0) {
    printConfigReference(std::cout, ConfigReferenceType::Markdown);
    return 0;
  }
  if (action == "xml" && args == 0) {
    printConfigReference(std::cout, ConfigReferenceType::XML);
    return 0;
  }
  if (action == "version" && args == 0) {
    std::cout << precice::getVersionInformation() << '\n';
    return 0;
  }

  printUsage();
  return 1;
}

#include <iostream>
#include <precice/Tooling.hpp>
#include <string>

void printUsage()
{
  std::cerr << "Usage:\n\n";
  std::cerr << "Print XML reference      :  precice-config-doc xml\n";
  std::cerr << "Print DTD for XML config :  precice-config-doc dtd\n";
  std::cerr << "Print Markdown reference :  precice-config-doc md\n";
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
    printConfigReference(std::cout, ConfigReferenceType::MD);
    return 0;
  }
  if (action == "xml" && args == 0) {
    printConfigReference(std::cout, ConfigReferenceType::XML);
    return 0;
  }

  printUsage();
  return 1;
}

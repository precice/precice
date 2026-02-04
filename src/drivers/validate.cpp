#include <iostream>
#include <precice/Exceptions.hpp>
#include <precice/Tooling.hpp>
#include <stdexcept>
#include <string>

void printUsage()
{
  std::cerr << "Usage:\n\n";
  std::cerr << "precice-config-validate FILE [ PARTICIPANT [ COMMSIZE ] ]\n";
}

int main(int argc, char **argv)
{
  const int args = argc - 1;

  if (args < 1 || args > 3) {
    printUsage();
    return 1;
  }

  using namespace precice::tooling;

  std::string file(argv[1]);
  std::string participant = (args > 1) ? std::string(argv[2]) : "";

  int size = 1;
  if (args == 3) {
    try {
      size = std::stoi(argv[3]);
    } catch (std::invalid_argument &e) {
      std::cerr << "ERROR: passed COMMSIZE is not a valid number\n";
      printUsage();
      return 1;
    }
  }

  try {
    precice::tooling::checkConfiguration(file, participant, size);
  } catch (const ::precice::Error &) {
    return 2;
  }
  return 0;
}

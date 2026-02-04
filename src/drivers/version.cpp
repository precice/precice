#include <iostream>
#include <precice/Tooling.hpp>

int main(int argc, char **argv)
{
  std::cout << precice::getVersionInformation() << '\n';
  return 0;
}

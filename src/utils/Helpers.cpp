#include "utils/Globals.hpp"

namespace precice {
namespace utils {

bool isMachineBigEndian()
{
   union {
      uint32_t i;
      char c[4];
   } bint = {0x01020304};

   return bint.c[0] == 1;
}

std::string getTypeName(const double& var)
{
  return std::string("float");
}

std::string getTypeName(const std::string& var)
{
  return std::string("string");
}

std::string getTypeName(const bool& var)
{
  return std::string("boolean");
}

std::string getTypeName(const int& var)
{
  return std::string("integer");
}

//template<>
//double getZero ( double ) { return 0.0; }
//
//template<>
//int getZero ( int ) { return 0; }

}}

#include "utils/Globals.hpp"

namespace precice {
namespace utils {


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

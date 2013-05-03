// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "utils/Globals.hpp"

namespace precice {
namespace utils {

void initializeZero ( double & toInitialize )
{
   toInitialize = 0.0;
}

void initializeZero ( int & toInitialize )
{
   toInitialize = 0;
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

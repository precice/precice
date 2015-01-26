// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "StringTest.hpp"
#include "../String.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include <string>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::StringTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log StringTest:: _log("precice::utils::StringTest");

StringTest:: StringTest()
:
  TestCase("utils::StringTest")
{}

void StringTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testWrapText);
    testMethod(testCheckAppendExtension);
  }
}

void StringTest:: testWrapText()
{
  preciceTrace("testWrapText()");
  std::string text("123456 1234567 12345678");
  std::string wrapped = wrapText(text, 4, 0);
  validateEquals(wrapped, std::string("123456\n1234567\n12345678"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 5, 0);
  validateEquals(wrapped, std::string("1234\n1234\n1234"));

  text = "1234 123 5";
  wrapped = wrapText(text, 5, 0);
  validateEquals(wrapped, std::string("1234\n123 5"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 7, 2);
  validateEquals(wrapped, std::string("1234\n  1234\n  1234"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 5, 2);
  validateEquals(wrapped, std::string("1234\n  1234\n  1234"));

  text = "12345678 1234 1 1234";
  wrapped = wrapText(text, 8, 2);
  validateEquals(wrapped, std::string("12345678\n  1234 1\n  1234"));
}

void StringTest:: testCheckAppendExtension()
{
  preciceTrace("testCheckAppendextension()");
  std::string filename("somefile");
  std::string extension(".xyz");

  std::string result(filename);
  checkAppendExtension(result, extension);
  validateWithMessage(result.compare(filename + extension) == 0, result);

  result = filename + extension;
  checkAppendExtension(result, extension);
  validateWithMessage(result.compare(filename + extension) == 0, result);

  result = filename + extension + ".zyx";
  checkAppendExtension(result, extension);
  validateWithMessage(result.compare(filename + extension + ".zyx" + extension) == 0,
                      result);
}


}}} // namespace precice, utils, tests


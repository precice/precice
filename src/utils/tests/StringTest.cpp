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
    testMethod(testTokenize);
    testMethod(testWrapText);
    testMethod(testCheckAppendExtension);
  }
}

void StringTest:: testTokenize()
{
  preciceTrace("tokenize()");

  std::vector<std::string> tokens;

  std::string stringA("This-is-a-test");
  std::string delimiterA("-");
  tokens = tokenize ( stringA, delimiterA);
  validateEquals(tokens.size(), 4);
  validateEquals(tokens[0], std::string("This"));
  validateEquals(tokens[1], std::string("is"));
  validateEquals(tokens[2], std::string("a"));
  validateEquals(tokens[3], std::string("test"));

  std::string stringB("This::is::another::test");
  std::string delimiterB("::");
  tokens = tokenize(stringB, delimiterB);
  validateEquals(tokens.size(), 4);
  validateEquals(tokens[0], std::string("This"));
  validateEquals(tokens[1], std::string("is"));
  validateEquals(tokens[2], std::string("another"));
  validateEquals(tokens[3], std::string("test"));

  std::string stringC("::Test");
  std::string delimiterC("::");
  tokens = tokenize(stringC, delimiterC);
  validateEquals(tokens.size(), 1 );
  validateEquals(tokens[0], std::string("Test"));

  std::string stringD("Test::");
  std::string delimiterD("::");
  tokens = tokenize(stringD, delimiterD);
  validateEquals(tokens.size(), 1);
  validateEquals(tokens[0], std::string("Test"));

  std::string stringE("::This::is:another::test::");
  std::string delimiterE("::");
  tokens = tokenize ( stringE, delimiterE);
  validateEquals(tokens.size(), 3);
  validateEquals(tokens[0], std::string("This"));
  validateEquals(tokens[1], std::string("is:another"));
  validateEquals(tokens[2], std::string("test"));

  {
    std::string string(" This is another test ");
    std::string delimiter(" ");
    tokens = tokenize(string, delimiter);
    validateEquals(tokens.size(), 4 );
    validateEquals(tokens[0], std::string("This"));
    validateEquals(tokens[1], std::string("is"));
    validateEquals(tokens[2], std::string("another"));
    validateEquals(tokens[3], std::string("test"));
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


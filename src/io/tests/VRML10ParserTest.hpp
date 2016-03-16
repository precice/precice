// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_TESTS_VRML10ParserTest_HPP_
#define PRECICE_IO_TESTS_VRML10ParserTest_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace io {
namespace tests {

/**
 * @brief Provides tests for class VRML10Parser.
 */
class VRML10ParserTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   VRML10ParserTest();

   /**
    * @brief Destructor.
    */
   virtual ~VRML10ParserTest() {}

   /**
    * @brief Empty.
    */
   virtual void setUp() {}

   /**
    * @brief Runs all tests, triggered by TestCaseCollection.
    */
   virtual void run();

private:

   // @brief Logging device.
   static logging::Logger _log;

   /**
    *
    */
   void testParseCube();
};

}}} // namespace precice, io, tests

#endif //PRECICE_IO_TESTS_VRML10ParserTest_HPP_

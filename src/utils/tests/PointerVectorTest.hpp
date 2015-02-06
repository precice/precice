// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_
#define PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Provides tests for class PointerVector.
 */
class PointerVectorTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   PointerVectorTest ();

   /**
    * @brief Destructor.
    */
   virtual ~PointerVectorTest() {};

   /**
    * Setup for tests, empty.
    */
   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_ */

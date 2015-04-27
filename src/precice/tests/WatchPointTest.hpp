// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_TESTS_WATCHPOINTTEST_HPP_
#define PRECICE_TESTS_WATCHPOINTTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace tests {

class WatchPointTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   WatchPointTest();

   /**
    * @brief Destructor.
    */
   virtual ~WatchPointTest() {}

   /**
    * @brief Empty.
    */
   virtual void setUp() {}

   /**
    * @brief Runs all tests of this class.
    */
   virtual void run();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}} // namespace precice, tests

#endif /* PRECICE_TESTS_WATCHPOINTTEST_HPP_ */

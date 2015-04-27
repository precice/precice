// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_TESTS_COMMUNICATEMESHTEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATEMESHTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace com {
namespace tests {

/**
 * @brief Provides tests for class GeometryHub.
 */
class CommunicateMeshTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   CommunicateMeshTest();

   /**
    * @brief Destructor.
    */
   virtual ~CommunicateMeshTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   /**
    * @brief Test with two involved solvers.
    */
   void testTwoSolvers ();

//   void testHubThreeSolvers ();
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_COMMUNICATEMESHTEST_HPP_ */

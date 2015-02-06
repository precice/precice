// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_
#define PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace action {
namespace tests {


class ModifyCoordinatesActionTest : public tarch::tests::TestCase
{
public:

   ModifyCoordinatesActionTest ();

   virtual ~ModifyCoordinatesActionTest () {};

   virtual void setUp () {}

   virtual void run ();

private:

   static tarch::logging::Log _log;

   /**
    * @brief Tests the setting of discplacements as absolute displacements.
    */
//   void testSetAsDisplacement ();

   /**
    * @brief Tests the adding of displacements as relative displacements.
    */
   void testAddToCoordinates ();

   void testSubtractFromCoordinates ();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_ */

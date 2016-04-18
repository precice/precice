// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>

namespace precice {
namespace geometry {
namespace tests {

class CuboidTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   CuboidTest ();

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Runs cuboid tests.
    */
   virtual void run ();

private:

   static logging::Logger _log;

   void testCreation ();

   void testConfiguration ();

   void testSubIDs2D ();

   void testSubIDs3D ();
};

}}} // namespace precice, geometry, tests

#endif // PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_

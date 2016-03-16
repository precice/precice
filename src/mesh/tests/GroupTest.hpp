// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_TESTS_GROUPTEST_HPP_
#define PRECICE_MESH_TESTS_GROUPTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Group.
 */
class GroupTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   GroupTest();

   /**
    * @Destructor.
    */
   virtual ~GroupTest() {}

   virtual void setUp() {}

   /**
    * @brief Runs all tests.
    */
   virtual void run();

private:

   // @brief Logging device.
   static logging::Logger _log;
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_GROUPTEST_HPP_

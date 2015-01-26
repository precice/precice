// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_TESTS_EDGEXTEST_HPP_
#define PRECICE_MESH_TESTS_EDGEXTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Edge.
 */
class EdgeTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   EdgeTest ();

   /**
    * @brief Destructor.
    */
   virtual ~EdgeTest () {};

   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   void test ();
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_EDGEXTEST_HPP_

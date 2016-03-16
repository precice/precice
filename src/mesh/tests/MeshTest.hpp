// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_TESTS_MESHTEST_HPP_
#define PRECICE_MESH_TESTS_MESHTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Mesh.
 */
class MeshTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   MeshTest ();

   /**
    * @brief Destructor.
    */
   virtual ~MeshTest () {}

   /**
    * @brief Setup for tests, empty.
    */
   virtual void setUp () {}

   /**
    * @brief Calls all test methods.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static logging::Logger _log;

   /**
    * @brief Tests creation of a simple mesh with vertices, edges, and triangle.
    */
   void testBasicSetup ();

   /**
    * @brief Tests assigning geometry sub-IDs to the mesh elements.
    */
   void testSubIDs ();

   /**
    * @brief Tests method Mesh::computeState().
    */
   void testComputeState ();

   void testBoundingBoxCOG();

   /**
    * @brief Demonstrates the capabilities of class Mesh.
    */
   void testDemonstration ();

   /**
    * @brief Tests method Mesh::computeDistribution
    */
   void testDistribution();
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_MESHTEST_HPP_

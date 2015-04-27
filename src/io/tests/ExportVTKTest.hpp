// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_TESTS_EXPORTVTKTEST_HPP_
#define PRECICE_IO_TESTS_EXPORTVTKTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace io {
namespace tests {

/**
 * @brief Provides tests for class ExportVTK.
 */
class ExportVTKTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
   ExportVTKTest();

   /**
    * @brief Destructor, empty.
    */
   virtual ~ExportVTKTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp() {}

   /**
    * @brief Runs all test methods.
    */
   virtual void run();

private:

   static tarch::logging::Log _log;

   /**
    * @brief Exports a two-dimensional polygonal surface mesh.
    */
   void testExportPolygonalMesh();

   /**
    * @brief Exports a three-dimensional triangulated surface mesh.
    */
   void testExportTriangulatedMesh();

   /**
    * @brief Exports a three-dimensional quadrialteral surface mesh.
    */
   void testExportQuadMesh();
};

}}} // namespace precice, io, tests

#endif // PRECICE_IO_TESTS_EXPORTVTKTEST_HPP_

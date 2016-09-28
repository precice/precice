#ifndef PRECICE_IO_TESTS_EXPORTVTKTEST_HPP_
#define PRECICE_IO_TESTS_EXPORTVTKTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

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

   static logging::Logger _log;

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

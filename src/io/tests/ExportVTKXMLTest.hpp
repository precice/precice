#ifndef PRECICE_IO_TESTS_EXPORTVTKXMLTEST_HPP_
#define PRECICE_IO_TESTS_EXPORTVTKXMLTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace io {
namespace tests {

/**
 * @brief Provides tests for class ExportVTKXML. The resemble closely the ExportVTK tests, but in parallel.
 */
class ExportVTKXMLTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
   ExportVTKXMLTest();

   /**
    * @brief Destructor, empty.
    */
   virtual ~ExportVTKXMLTest() {};

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

   void setUpMasterSlave();

   void tearDownMasterSlave();
};

}}} // namespace precice, io, tests

#endif // PRECICE_IO_TESTS_EXPORTVTKXMLTEST_HPP_

#ifndef PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace geometry {
namespace tests {

/**
 * @brief Provides tests for class ImportGeometry.
 */
class ImportGeometryTest : public tarch::tests::TestCase
{
public:

   ImportGeometryTest ();

   virtual ~ImportGeometryTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static tarch::logging::Log _log;

   void testImportVRMLConfig ();
};

}}} // namespace precice, geometry, tests

#endif /* PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_ */

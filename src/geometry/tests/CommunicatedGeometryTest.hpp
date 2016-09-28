#ifndef PRECICE_NO_MPI
 
#ifndef PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace geometry {
namespace tests {

/**
 * @brief Provides tests for class CommunicatedGeometry.
 */
class CommunicatedGeometryTest : public tarch::tests::TestCase
{
public:

  CommunicatedGeometryTest ();

   virtual ~CommunicatedGeometryTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;

   void testScatterMesh ();

   void testGatherMesh ();
};

}}} // namespace precice, geometry, tests

#endif /* PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_ */
#endif // PRECICE_NO_MPI

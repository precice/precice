#ifndef PRECICE_NO_MPI
 
#ifndef PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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

   static tarch::logging::Log _log;

   void testScatterMesh ();

   void testGatherMesh ();
};

}}} // namespace precice, geometry, tests

#endif /* PRECICE_GEOMETRY_TESTS_COMMNUICATEDGEOMETRYTEST_HPP_ */
#endif // PRECICE_NO_MPI

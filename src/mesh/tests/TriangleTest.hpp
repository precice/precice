#ifndef PRECICE_MESH_TESTS_TRIANGLETEST_HPP_
#define PRECICE_MESH_TESTS_TRIANGLETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Triangle.
 */
class TriangleTest
:
   public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   TriangleTest ();

   /**
    * @brief Destructor.
    */
   virtual ~TriangleTest () {};

   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   /**
    * @brief Test.
    */
   void test ();
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_TRIANGLETEST_HPP_

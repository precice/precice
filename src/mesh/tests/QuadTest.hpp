#ifndef PRECICE_MESH_TESTS_QUADTEST_HPP_
#define PRECICE_MESH_TESTS_QUADTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Quad.
 */
class QuadTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
  QuadTest ();

   /**
    * @brief Destructor.
    */
   virtual ~QuadTest () {};

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

#endif // PRECICE_MESH_TESTS_QUADTEST_HPP_

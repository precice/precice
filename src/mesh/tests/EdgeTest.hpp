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

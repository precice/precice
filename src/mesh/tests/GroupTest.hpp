#ifndef PRECICE_MESH_TESTS_GROUPTEST_HPP_
#define PRECICE_MESH_TESTS_GROUPTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Group.
 */
class GroupTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   GroupTest();

   /**
    * @Destructor.
    */
   virtual ~GroupTest() {}

   virtual void setUp() {}

   /**
    * @brief Runs all tests.
    */
   virtual void run();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_GROUPTEST_HPP_

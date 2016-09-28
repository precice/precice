#ifndef PRECICE_MESH_TESTS_PROPERTYMESHTEST_HPP_
#define PRECICE_MESH_TESTS_PROPERTYMESHTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
namespace tests {


/**
 * @brief Provides tests for class PropertyContainer.
 */
class PropertyContainerTest
:
   public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   PropertyContainerTest ();

   /*
    * @brief Destructor.
    */
   virtual ~PropertyContainerTest() {};

   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static logging::Logger _log;

   /**
    * @brief Tests a single property.
    */
   void testSinglePropertyContainer ();

   /**
    * @brief Tests two properties connected in a hierarchical relationship.
    */
   void testHierarchicalPropertyContainers ();

   void testMultipleParents ();
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_PROPERTYMESHTEST_HPP_

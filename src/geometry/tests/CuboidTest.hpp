#ifndef PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>

namespace precice {
namespace geometry {
namespace tests {

class CuboidTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   CuboidTest ();

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Runs cuboid tests.
    */
   virtual void run ();

private:

   static logging::Logger _log;

   void testCreation ();

   void testConfiguration ();

   void testSubIDs2D ();

   void testSubIDs3D ();
};

}}} // namespace precice, geometry, tests

#endif // PRECICE_GEOMETRY_TESTS_CUBOIDTEST_HPP_

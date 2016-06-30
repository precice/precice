#ifndef PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_
#define PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace action {
namespace tests {


class ModifyCoordinatesActionTest : public tarch::tests::TestCase
{
public:

   ModifyCoordinatesActionTest ();

   virtual ~ModifyCoordinatesActionTest () {};

   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;

   /**
    * @brief Tests the setting of discplacements as absolute displacements.
    */
//   void testSetAsDisplacement ();

   /**
    * @brief Tests the adding of displacements as relative displacements.
    */
   void testAddToCoordinates ();

   void testSubtractFromCoordinates ();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_TESTS_MODIFYCOORDINATESACTIONTEST_HPP_ */

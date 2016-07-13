#ifndef PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_
#define PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Provides tests for class PointerVector.
 */
class PointerVectorTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   PointerVectorTest ();

   /**
    * @brief Destructor.
    */
   virtual ~PointerVectorTest() {};

   /**
    * Setup for tests, empty.
    */
   virtual void setUp () {}

   /**
    * @brief Runs all tests.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_POINTERVECTORTEST_HPP_ */

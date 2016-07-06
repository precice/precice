#ifndef PRECICE_TESTS_WATCHPOINTTEST_HPP_
#define PRECICE_TESTS_WATCHPOINTTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace tests {

class WatchPointTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   WatchPointTest();

   /**
    * @brief Destructor.
    */
   virtual ~WatchPointTest() {}

   /**
    * @brief Empty.
    */
   virtual void setUp() {}

   /**
    * @brief Runs all tests of this class.
    */
   virtual void run();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}} // namespace precice, tests

#endif /* PRECICE_TESTS_WATCHPOINTTEST_HPP_ */

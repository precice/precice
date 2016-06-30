#ifndef PRECICE_UTILS_TESTS_HELPERSTEST_HPP_
#define PRECICE_UTILS_TESTS_HELPERSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Tests methods from file Helpers.
 */
class HelpersTest : public tarch::tests::TestCase
{
public:

   HelpersTest ();

   virtual void setUp () {}

   virtual void run ();

private:

    static logging::Logger _log;

//    void testGetZero ();

    void testAppendTo ();

    void testOperatorPlusForVectors ();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_HELPERSTEST_HPP_ */

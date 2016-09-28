#ifndef PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_
#define PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {
namespace tests {


class ManageUniqueIDsTest : public tarch::tests::TestCase
{
public:

  ManageUniqueIDsTest();

  virtual ~ManageUniqueIDsTest() {};

  virtual void setUp () {}

  virtual void run ();

private:

  static logging::Logger _log;

  void testUniqueIDs ();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_ */

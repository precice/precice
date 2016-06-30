#ifndef PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_
#define PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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

  static tarch::logging::Log _log;

  void testUniqueIDs ();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_MANAGEUNIQUEIDSTEST_HPP_ */

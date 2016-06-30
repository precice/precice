#ifndef PRECICE_UTILS_DIMENSIONSTEST_HPP_
#define PRECICE_UTILS_DIMENSIONSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace utils {
namespace tests {

/**
 * Provides tests for methods/classes in file utils/Dimensions.
 */
class DimensionsTest : public tarch::tests::TestCase
{
public:

  DimensionsTest();

  virtual ~DimensionsTest() {}

  virtual void setUp() {}

  virtual void run();

private:

  static tarch::logging::Log _log;

  void testLinearizeDelinearize();

//  void testGetHyperfaceCornerIndices();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_DIMENSIONSTEST_HPP_ */

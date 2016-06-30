#ifndef PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_
#define PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mapping {
namespace tests {

/**
 * @brief Provides tests for class NearestNeighborMapping.
 */
class NearestNeighborMappingTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  NearestNeighborMappingTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~NearestNeighborMappingTest() {}

  /**
   * @brief Prepares run of tests, empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  void testConsistentNonIncremental();

  void testConservativeNonIncremental();
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_ */

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_M2N_TESTS_POINT_TO_POINT_COMMUNICATION_TEST_HPP_
#define PRECICE_M2N_TESTS_POINT_TO_POINT_COMMUNICATION_TEST_HPP_

#include "com/CommunicationFactory.hpp"

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace m2n {
namespace tests {

/**
 * @brief Provides tests for class PointToPointCommunication.
 */
class PointToPointCommunicationTest : public tarch::tests::TestCase {
public:
  PointToPointCommunicationTest();

  virtual ~PointToPointCommunicationTest(){};

  /**
   * @brief Empty.
   */
  virtual void
  setUp() {
  }

  virtual void run();

private:
  static tarch::logging::Log _log;

  void testSocketCommunication();

  void testMPIPortsCommunication();

  void test(com::CommunicationFactory::SharedPointer cf);
};
}
}
} // namespace precice, m2n, tests

#endif /* PRECICE_M2N_TESTS_POINT_TO_POINT_COMMUNICATION_TEST_HPP_ */

#endif // not PRECICE_NO_MPI

#ifndef PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace com {
namespace tests {

class MPIDirectCommunicationTest : public tarch::tests::TestCase
{
public:

  MPIDirectCommunicationTest();

  virtual ~MPIDirectCommunicationTest() {}

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  virtual void run();

private:

  static logging::Logger _log;

# ifndef PRECICE_NO_MPI

  void testSendReceiveTwoProcesses();

  void testSendReceiveThreeProcesses();

# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_ */

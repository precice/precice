#ifndef PRECICE_NO_MPI

#ifndef PRECICE_M2N_TESTS_GATHER_SCATTER_COMMUNICATION_TEST_HPP_
#define PRECICE_M2N_TESTS_GATHER_SCATTER_COMMUNICATION_TEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace m2n {
namespace tests {

/**
 * @brief Provides tests for class GatherScatterCommunication.
 */
class GatherScatterCommunicationTest : public tarch::tests::TestCase
{
public:

  GatherScatterCommunicationTest ();

   virtual ~GatherScatterCommunicationTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;

   /**
    * @brief Tests the receive and send of a double array from a serial participant
    * to a participant running on 3 processes
    */
   void testSendReceiveAll ();
};

}}} // namespace precice, m2n, tests

#endif /* PRECICE_M2N_TESTS_GATHER_SCATTER_COMMUNICATION_TEST_HPP_ */
#endif // PRECICE_NO_MPI

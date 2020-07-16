#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIPorts,
                      *boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(SendAndReceiveMM)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::mastermaster;
  TestSendAndReceive<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendAndReceiveMS)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::masterslave;
  TestSendAndReceive<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesMM)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::mastermaster;
  TestSendReceiveFourProcesses<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // MPIPortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

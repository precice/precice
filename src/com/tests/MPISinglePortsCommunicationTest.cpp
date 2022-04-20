#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPISinglePortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPISinglePorts,
                      *boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(SendAndReceiveMM)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceive<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendAndReceiveMS)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::primarysecondary;
  TestSendAndReceive<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesMM)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendReceiveFourProcesses<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // MPISinglePortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

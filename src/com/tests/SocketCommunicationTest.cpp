#include <vector>
#include "GenericTestFunctions.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketCommunication.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_TEST_SPECIALIZED_COLLECTION_COMPARE(std::vector<int>)

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(Socket)

BOOST_AUTO_TEST_CASE(SendAndReceiveMM)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::mastermaster;
  TestSendAndReceive<SocketCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendAndReceiveMS)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::masterslave;
  TestSendAndReceive<SocketCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::mastermaster;
  TestSendReceiveFourProcesses<SocketCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<SocketCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<SocketCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<SocketCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

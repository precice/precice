#include "GenericTestFunctions.hpp"
#include "com/SocketCommunication.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_TEST_SPECIALIZED_COLLECTION_COMPARE(std::vector<int>)

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(Socket)

BOOST_AUTO_TEST_CASE(SendAndReceive)

{
  PRECICE_TEST(2_ranks, Require::Events);
  TestSendAndReceive<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcesses<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient)
{
  PRECICE_TEST(2_ranks, Require::Events);
  TestSendReceiveTwoProcessesServerClient<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient)
{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcessesServerClient<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2)
{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcessesServerClientV2<SocketCommunication>();
}

BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

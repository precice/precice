#include "com/SocketCommunication.hpp"
#include "testing/Testing.hpp"
#include "GenericTestFunctions.hpp"

using namespace precice;
using namespace precice::com;


BOOST_TEST_SPECIALIZED_COLLECTION_COMPARE(std::vector<int>)

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(Socket)


BOOST_AUTO_TEST_CASE(SendAndReceive,
                     * testing::MinRanks(2))
{
  TestSendAndReceive<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
{
  TestSendReceiveFourProcesses<SocketCommunication>();
}

/// Disabled because acceptConnection lost ability to accept multiple connections
BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesV2,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::disabled())
{
  TestSendReceiveFourProcessesV2<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
{
  TestSendReceiveTwoProcessesServerClient<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
{
  TestSendReceiveFourProcessesServerClient<SocketCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
{
  TestSendReceiveFourProcessesServerClientV2<SocketCommunication>();
}


BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

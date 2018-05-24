#include "com/SocketCommunication.hpp"
#include "testing/Testing.hpp"
#include "SendAndReceive.hpp"

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

BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

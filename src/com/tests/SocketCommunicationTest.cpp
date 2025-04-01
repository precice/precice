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

BOOST_AUTO_TEST_SUITE(Intra)

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceivePrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveRanges)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveRanges<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveEigen<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastPrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastVectors<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReducePrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReduceVectors<SocketCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Intra

BOOST_AUTO_TEST_SUITE(Inter)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceivePrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveEigen<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveRanges)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveRanges<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastPrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastVectors<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestReducePrimitiveTypes<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestReduceVectors<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendReceiveFourProcesses<SocketCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Inter

BOOST_AUTO_TEST_SUITE(Server)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveTwo)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFour)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<SocketCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFourV2)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<SocketCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Server

BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

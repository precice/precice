#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPISinglePortsCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPISinglePorts)

BOOST_AUTO_TEST_SUITE(Intra)

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceivePrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveEigen<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastPrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastVectors<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReducePrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReduceVectors<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Intra

BOOST_AUTO_TEST_SUITE(Inter)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceivePrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveRanges)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveRanges<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveEigen<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastPrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastEigen<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestReducePrimitiveTypes<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestReduceVectors<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST();
  using namespace precice::testing::com::primaryprimary;
  TestSendReceiveFourProcesses<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Inter

BOOST_AUTO_TEST_SUITE(Server)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveTwo)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFour)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<MPISinglePortsCommunication>(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveFourV2)
{
  PRECICE_TEST();
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<MPISinglePortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Server

BOOST_AUTO_TEST_SUITE_END() // MPISinglePorts
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif

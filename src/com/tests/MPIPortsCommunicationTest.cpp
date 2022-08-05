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

BOOST_AUTO_TEST_SUITE(MPIPorts)

BOOST_AUTO_TEST_SUITE(Intra)

BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestSendAndReceivePrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveRanges)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveRanges<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveEigen<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestBroadcastPrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestBroadcastVectors<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestReducePrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestReduceVectors<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Intra

BOOST_AUTO_TEST_SUITE(Inter)

BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceivePrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveRanges)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveRanges<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendAndReceiveEigen<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastPrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestBroadcastVectors<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestReducePrimitiveTypes<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestReduceVectors<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::primaryprimary;
  TestSendReceiveFourProcesses<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Inter

BOOST_AUTO_TEST_SUITE(Server)

BOOST_AUTO_TEST_CASE(SendReceiveTwo)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveTwoProcessesServerClient<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFour)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClient<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveFourV2)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks), Require::Events);
  using namespace precice::testing::com::serverclient;
  TestSendReceiveFourProcessesServerClientV2<MPIPortsCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Server

BOOST_AUTO_TEST_SUITE_END() // MPIPorts
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif

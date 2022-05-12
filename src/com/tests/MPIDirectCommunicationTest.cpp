#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIDirect)

BOOST_AUTO_TEST_SUITE(Intra)

BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestSendAndReceivePrimitiveTypes<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_CASE(SendReceiveVectors)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveVectors<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestBroadcastPrimitiveTypes<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_CASE(BroadcastVectors)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestBroadcastVectors<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestReducePrimitiveTypes<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_CASE(ReduceVectors)
{
  PRECICE_TEST(2_ranks, Require::Events);
  using namespace precice::testing::com::intracomm;
  TestReduceVectors<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Intra

BOOST_AUTO_TEST_SUITE_END() // MPIDirect
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif

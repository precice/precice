#ifndef PRECICE_NO_MPI

#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "com/tests/GenericTestFunctions.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using precice::com::MPIDirectCommunication;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIDirect)

BOOST_AUTO_TEST_SUITE(Intra)

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceivePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceivePrimitiveTypes<MPIDirectCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestSendAndReceiveEigen<MPIDirectCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastPrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastPrimitiveTypes<MPIDirectCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestBroadcastVectors<MPIDirectCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReducePrimitives)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReducePrimitiveTypes<MPIDirectCommunication>(context);
}

PRECICE_TEST_SETUP(2_ranks, Require::Events)
BOOST_AUTO_TEST_CASE(ReduceEigen)
{
  PRECICE_TEST();
  using namespace precice::testing::com::intracomm;
  TestReduceVectors<MPIDirectCommunication>(context);
}

BOOST_AUTO_TEST_SUITE_END() // Intra

BOOST_AUTO_TEST_SUITE_END() // MPIDirect
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif

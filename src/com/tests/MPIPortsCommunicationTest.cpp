#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIPorts,
                      *boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(SendAndReceive)
{
  PRECICE_TEST(2_ranks, Require::Events);
  TestSendAndReceive<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcesses)
{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcesses<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient)

{
  PRECICE_TEST(2_ranks, Require::Events);
  TestSendReceiveTwoProcessesServerClient<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient)

{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcessesServerClient<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2)
{
  PRECICE_TEST(4_ranks, Require::Events);
  TestSendReceiveFourProcessesServerClientV2<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_SUITE_END() // MPIPortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

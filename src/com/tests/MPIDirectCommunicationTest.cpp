#ifndef PRECICE_NO_MPI

#include "GenericTestFunctions.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIDirect)

BOOST_AUTO_TEST_CASE(SendAndReceive,
                     *testing::MinRanks(2) * boost::unit_test::fixture<testing::SyncProcessesFixture>() * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  TestSendAndReceive<MPIDirectCommunication>();
}

BOOST_AUTO_TEST_SUITE_END() // MPIDirectCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

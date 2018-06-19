#ifndef PRECICE_NO_MPI

#include "com/MPIPortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"
#include "GenericTestFunctions.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIPorts,
                      * boost::unit_test::label("MPI_Ports"))


BOOST_AUTO_TEST_CASE(SendAndReceive,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  TestSendAndReceive<MPIPortsCommunication>();
}


BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
       
{
  TestSendReceiveTwoProcessesServerClient<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
       
{
  TestSendReceiveFourProcessesServerClient<MPIPortsCommunication>();
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClientV2,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
{
  TestSendReceiveFourProcessesServerClientV2<MPIPortsCommunication>();
}


BOOST_AUTO_TEST_SUITE_END() // MPIPortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "com/MPIDirectCommunication.hpp"
#include "testing/Testing.hpp"
#include "GenericTestFunctions.hpp"


using Par = precice::utils::Parallel;
using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIDirect,
                      * testing::MinRanks(2)
                      * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                      * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))

BOOST_AUTO_TEST_CASE(SendAndReceive)
{
  if (Par::getCommunicatorSize() != 2)
    return;

  TestSendAndReceive<MPIDirectCommunication>();
}

BOOST_AUTO_TEST_SUITE_END() // MPIDirectCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

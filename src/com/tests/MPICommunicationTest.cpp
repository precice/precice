#ifndef PRECICE_NO_MPI

#include "com/MPIPortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "SendAndReceive.hpp"

using Par = precice::utils::Parallel;
using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPICommunication,
                      * testing::MinRanks(2)
                      * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1}))
                      * boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(SendAndReceive)
{
  TestSendAndReceive<MPIPortsCommunication>();
}


BOOST_AUTO_TEST_SUITE_END() // MPICommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "../MPIDirectCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using Par = precice::utils::Parallel;
using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIDirect)

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcesses,
                     *testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1}))
                     * boost::unit_test::label("MPI_Ports"))
{
  if (Par::getCommunicatorSize() != 2)
    return;

  if (Par::getProcessRank() < 2) {
    MPIDirectCommunication communication;
    std::string            nameEven("even");
    std::string            nameOdd("odd");

    if (Par::getProcessRank() == 0) {
      Par::splitCommunicator(nameEven);
      communication.acceptConnection(nameEven, nameOdd);
      int message = 1;
      communication.send(message, 0);
      communication.receive(message, 0);
      BOOST_TEST(message == 2);
      communication.closeConnection();
    } else if (Par::getProcessRank() == 1) {
      Par::splitCommunicator(nameOdd);
      communication.requestConnection(nameEven, nameOdd, 0, 1);
      int message = -1;
      communication.receive(message, 0);
      BOOST_TEST(message == 1);
      message = 2;
      communication.send(message, 0);
      communication.closeConnection();
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // MPIDirectCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

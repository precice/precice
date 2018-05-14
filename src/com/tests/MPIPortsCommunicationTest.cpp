#ifndef PRECICE_NO_MPI

#include "com/MPIPortsCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MPIPorts,
                      * boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcesses,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  MPIPortsCommunication communication;

  std::string nameEven("even");
  std::string nameOdd("odd");

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    communication.acceptConnection(nameEven, nameOdd);
    int message = 1;
    communication.send(message, 0);
    communication.receive(message, 0);
    BOOST_TEST(message == 2);
    communication.closeConnection();
    break;
  }
  case 1: {
    communication.requestConnection(nameEven, nameOdd, 0, 1);
    int message = -1;
    communication.receive(message, 0);
    BOOST_TEST(message == 1);
    message = 2;
    communication.send(message, 0);
    communication.closeConnection();
    break;
  }
  }
}

BOOST_AUTO_TEST_CASE(SendReceiveTwoProcessesServerClient,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>()
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
       
{
  MPIPortsCommunication communication;

  std::string nameEven("even");
  std::string nameOdd("odd");

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    communication.acceptConnectionAsServer(nameEven, nameOdd, 1);
    int message = 1;
    communication.send(message, 0);
    communication.receive(message, 0);
    BOOST_TEST(message == 2);
    communication.closeConnection();
    break;
  }
  case 1: {
    communication.requestConnectionAsClient(nameEven, nameOdd);
    int message = -1;
    communication.receive(message, 0);
    BOOST_TEST(message == 1);
    message = 2;
    communication.send(message, 0);
    communication.closeConnection();
    break;
  }
  }
}

BOOST_AUTO_TEST_CASE(SendReceiveFourProcessesServerClient,
                     * testing::MinRanks(4)
                     * boost::unit_test::fixture<testing::SyncProcessesFixture>())
       
{
  MPIPortsCommunication communication;

  std::string nameEven("even");
  std::string nameOdd("odd");

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    communication.acceptConnectionAsServer(nameEven, nameOdd, 2);

    int requestorLocalRank = 0;
    int requestorGlobalRank = -1;
    communication.receive(requestorGlobalRank, requestorLocalRank);
    BOOST_TEST(requestorGlobalRank >= 2);
    BOOST_TEST(requestorGlobalRank <= 3);
    int message = requestorGlobalRank * 10;
    communication.send(message, requestorLocalRank);
    communication.receive(message, requestorLocalRank);
    BOOST_TEST(message == requestorGlobalRank * 10 + 2);

    requestorLocalRank = 1;
    requestorGlobalRank = -1;
    communication.receive(requestorGlobalRank, requestorLocalRank);
    BOOST_TEST(requestorGlobalRank >= 2);
    BOOST_TEST(requestorGlobalRank <= 3);
    message = requestorGlobalRank * 10;
    communication.send(message, requestorLocalRank);
    communication.receive(message, requestorLocalRank);
    BOOST_TEST(message == requestorGlobalRank * 10 + 2);

    communication.closeConnection();
    break;
  }
  case 1: {
    // does not accept a connection
    break;
  }
  case 2: {
    communication.requestConnectionAsClient(nameEven, nameOdd);
    int globalRank = 2;
    communication.send(globalRank,0);

    int message = -1;
    communication.receive(message, 0);
    BOOST_TEST(message == 20);
    message += 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  case 3: {
    communication.requestConnectionAsClient(nameEven, nameOdd);
    int globalRank = 3;
    communication.send(globalRank,0);

    int message = -1;
    communication.receive(message, 0);
    BOOST_TEST(message == 30);
    message += 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  }
}


BOOST_AUTO_TEST_SUITE_END() // MPIPortsCommunication

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI

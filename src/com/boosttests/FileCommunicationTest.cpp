#include "com/FileCommunication.hpp"
#include "math/math.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(File,
                      *testing::OnRanks({0, 1}) * boost::unit_test::disabled())

BOOST_AUTO_TEST_CASE(SimpleSendReceive)
{
  int rank = utils::Parallel::getProcessRank();
  bool binaryMode = false;
  FileCommunication comTxt(binaryMode, "");
  BOOST_TEST(comTxt.isConnected() == false);
  std::string requester("FileCommunicationTest-testSimpleSendReceive-Requester-txt");
  std::string acceptor("FileCommunicationTest-testSimpleSendReceive-Acceptor-txt");
  Eigen::Vector3d doubleVector = Eigen::Vector3d::Constant(0.1234567890123456);
  Eigen::Matrix<int, 5, 1> intVector;
  intVector << 1, 2, 3, 4, 5;
  // int intVector[] = { 1, 2, 3, 4, 5 };
  int integer = 1;
  double dbl = 1.0;
  bool boolean = true;
  if (rank == 0) {
    comTxt.requestConnection(requester, acceptor, 0, 1);
    BOOST_CHECK(comTxt.isConnected());
    comTxt.startSendPackage(0);
    comTxt.send(std::string("rank0"), 0);
    comTxt.send(doubleVector.data(), doubleVector.size(), 0);
    comTxt.send(intVector.data(), 5, 0);
    comTxt.send(dbl, 0);
    comTxt.send(integer, 0);
    comTxt.send(boolean, 0);
    comTxt.finishSendPackage();
  } else {
    assertion(rank == 1, rank);
    comTxt.acceptConnection(requester, acceptor, 0, 1);
    BOOST_CHECK(comTxt.isConnected());
    comTxt.startReceivePackage(0);
    std::string message("");
    comTxt.receive(message, 0);
    double rawVec[3];
    Eigen::Map<Eigen::Vector3d>(rawVec).setConstant(0);
    comTxt.receive(rawVec, 3, 0);
    int intArrayReceived[] = {0, 0, 0, 0, 0};
    comTxt.receive(intArrayReceived, 5, 0);
    double dblReceived = 0.0;
    comTxt.receive(dblReceived, 0);
    int integerReceived = 0;
    comTxt.receive(integerReceived, 0);
    bool booleanReceived = false;
    comTxt.receive(booleanReceived, 0);
    comTxt.finishReceivePackage();
    BOOST_TEST(message == std::string("rank0"));
    BOOST_TEST(testing::equals(Eigen::Map<Eigen::Vector3d>(rawVec), doubleVector));
    BOOST_TEST(testing::equals(Eigen::Map<Eigen::Matrix<int, 5, 1>>(intArrayReceived), intVector));
    BOOST_TEST(testing::equals(dblReceived, dbl));
    BOOST_TEST(integerReceived == integer);
    BOOST_TEST(booleanReceived == boolean);
  }
  comTxt.closeConnection();
  BOOST_CHECK(not comTxt.isConnected());

  binaryMode = true;
  FileCommunication comBin(binaryMode, "");
  BOOST_TEST(comBin.isConnected() == false);
  requester = "FileCommunicationTest-testSimpleSendReceive-Requester-bin";
  acceptor = "FileCommunicationTest-testSimpleSendReceive-Acceptor-bin";
  if (rank == 0) {
    comBin.requestConnection(requester, acceptor, 0, 1);
    BOOST_CHECK(comBin.isConnected());
    comBin.startSendPackage(0);
    comBin.send(std::string("rank0"), 0);
    comBin.send(doubleVector.data(), 3, 0);
    comBin.send(intVector.data(), 5, 0);
    comBin.send(dbl, 0);
    comBin.send(integer, 0);
    comBin.send(boolean, 0);
    comBin.finishSendPackage();
  } else {
    assertion(rank == 1, rank);
    comBin.acceptConnection(requester, acceptor, 0, 1);
    BOOST_CHECK(comBin.isConnected());
    comBin.startReceivePackage(0);
    std::string message("");
    comBin.receive(message, 0);
    double rawVec[3];
    Eigen::Map<Eigen::Vector3d>(rawVec).setConstant(0);
    comBin.receive(rawVec, 3, 0);
    int intArrayReceived[] = {0, 0, 0, 0, 0};
    comBin.receive(intArrayReceived, 5, 0);
    double dblReceived = 0.0;
    comBin.receive(dblReceived, 0);
    int integerReceived = 0;
    comBin.receive(integerReceived, 0);
    bool booleanReceived = false;
    comBin.receive(booleanReceived, 0);
    comBin.finishReceivePackage();
    BOOST_TEST(message == std::string("rank0"));
    BOOST_TEST(testing::equals(Eigen::Map<Eigen::Vector3d>(rawVec), doubleVector));
    BOOST_TEST(testing::equals(Eigen::Map<Eigen::Matrix<int, 5, 1>>(intArrayReceived), intVector));
    BOOST_TEST(testing::equals(dblReceived, dbl));
    BOOST_TEST(integerReceived, integer);
    BOOST_TEST(booleanReceived, boolean);
  }
  comBin.closeConnection();
  BOOST_CHECK(not comBin.isConnected());
}

BOOST_AUTO_TEST_CASE(testMultipleExchanges)
{
  int rank = utils::Parallel::getProcessRank();
  bool binaryMode = false;
  FileCommunication com(binaryMode, "");
  BOOST_CHECK(not com.isConnected());
  std::string requester("FileCommunicationTest-testMultipleExchanges-Requester-txt");
  std::string acceptor("FileCommunicationTest-testMultipleExchanges-Acceptor-txt");
  int value0 = 1;
  int value1 = 2;
  int value2 = 3;
  if (rank == 0) {
    com.requestConnection(acceptor, requester, 0, 1);
    BOOST_CHECK(com.isConnected());

    com.startSendPackage(0);
    com.send(value0, 0);
    com.finishSendPackage();

    com.startReceivePackage(0);
    int number = 0;
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value1);

    com.startSendPackage(0);
    com.send(value2, 0);
    com.finishSendPackage();

    com.startSendPackage(0);
    com.send(value0, 0);
    com.finishSendPackage();

    com.startSendPackage(0);
    com.send(value1, 0);
    com.finishSendPackage();

    com.startReceivePackage(0);
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value2);

    com.closeConnection();
    BOOST_CHECK(not com.isConnected());
  } else {
    assertion(rank == 1, rank);
    com.acceptConnection(acceptor, requester, 0, 1);
    BOOST_CHECK(com.isConnected());

    com.startReceivePackage(0);
    int number = 0;
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value0);

    com.startSendPackage(0);
    com.send(value1, 0);
    com.finishSendPackage();

    com.startReceivePackage(0);
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value2);

    com.startReceivePackage(0);
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value0);

    com.startReceivePackage(0);
    com.receive(number, 0);
    com.finishReceivePackage();
    BOOST_TEST(number == value1);

    com.startSendPackage(0);
    com.send(value2, 0);
    com.finishSendPackage();

    com.closeConnection();
    BOOST_CHECK(not com.isConnected());
  }
}

BOOST_AUTO_TEST_SUITE_END() // FileCommunication
BOOST_AUTO_TEST_SUITE_END() // Communication

#include "com/SocketCommunication.hpp"
#include "testing/Testing.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(Socket)

BOOST_AUTO_TEST_CASE(SendAndReceive,
                     *testing::OnSize(2))
{
  SocketCommunication com;
  if (utils::Parallel::getProcessRank() == 0) {
    com.acceptConnection("process0", "process1");
    {
      std::string msg("testOne");
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_TEST(msg == std::string("testTwo"));
    }
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.receive(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(1)));
      msg = Eigen::Vector3d::Constant(2);
      com.send(msg.data(), msg.size(), 0);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.receive(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(1)));
      msg = Eigen::Vector4i::Constant(0);
      com.send(msg.data(), msg.size(), 0);
    }
    {
      std::vector<int> msg;
      com.receive(msg, 0);
      BOOST_CHECK(msg == std::vector<int>({1, 2, 3}));
      com.send(msg, 0);
    }
    {
      std::vector<double> msg;
      com.receive(msg, 0);
      BOOST_CHECK(msg == std::vector<double>({1.1, 2.2, 3.3}));
      com.send(msg, 0);
    }
    {
      double msg = 0.0;
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_TEST(msg == 1.0);
    }
    {
      int msg = 1;
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_TEST(msg == 2);
    }
    {
      bool msg = true;
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_TEST(msg == false);
    }
    com.closeConnection();
  } else if (utils::Parallel::getProcessRank() == 1) {
    com.requestConnection("process0", "process1", 0, 1);
    {
      std::string msg;
      com.receive(msg, 0);
      BOOST_TEST(msg == std::string("testOne"));
      msg = "testTwo";
      com.send(msg, 0);
    }
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(1);
      com.send(msg.data(), msg.size(), 0);
      com.receive(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(2)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(1);
      com.send(msg.data(), msg.size(), 0);
      com.receive(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Zero()));
    }
    {
      std::vector<int> msg{1, 2, 3};
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(msg == std::vector<int>({1, 2, 3}));
    }
    {
      std::vector<double> msg{1.1, 2.2, 3.3};
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(msg == std::vector<double>({1.1, 2.2, 3.3}));
    }
    {
      double msg = 1.0;
      com.receive(msg, 0);
      BOOST_TEST(msg == 0.0);
      msg = 1.0;
      com.send(msg, 0);
    }
    {
      int msg = 0;
      com.receive(msg, 0);
      BOOST_TEST(msg == 1);
      msg = 2;
      com.send(msg, 0);
    }
    {
      bool msg = false;
      com.receive(msg, 0);
      BOOST_TEST(msg == true);
      msg = false;
      com.send(msg, 0);
    }
    com.closeConnection();
  }
}

BOOST_AUTO_TEST_CASE(ParallelClient,
                     *testing::OnSize(3))
{
  SocketCommunication com;
  int                 rank = utils::Parallel::getProcessRank();
  if (rank == 0) {
    com.acceptConnection("server", "client");
    BOOST_TEST(com.getRemoteCommunicatorSize() == 2);
    std::string msg;
    com.receive(msg, 0);
    BOOST_TEST(msg == std::string("process 0"));
    com.receive(msg, 1);
    BOOST_TEST(msg == std::string("process 1"));
    int sendMsg = 1;
    com.send(sendMsg, 0);
    sendMsg = 2;
    com.send(sendMsg, 1);
    com.closeConnection();
  } else if ((rank == 1) || (rank == 2)) {
    com.requestConnection("server", "client", rank - 1, 2);
    BOOST_TEST(com.getRemoteCommunicatorSize() == 1);
    std::ostringstream rankMsg;
    rankMsg << "process " << rank - 1;
    com.send(rankMsg.str(), 0);
    int receiveMsg = 0;
    if (rank == 1) {
      com.receive(receiveMsg, 0);
      BOOST_TEST(receiveMsg == 1);
    } else {
      assertion(rank == 2, rank);
      com.receive(receiveMsg, 0);
      BOOST_TEST(receiveMsg == 2);
    }
    com.closeConnection();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Socket
BOOST_AUTO_TEST_SUITE_END() // Communication

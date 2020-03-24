#pragma once

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <string>
#include <vector>
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

/// Generic test function that is called from the tests for MPIPortsCommunication,
/// MPIDirectCommunication and SocketCommunication

namespace precice {
namespace testing {
namespace com {

namespace mastermaster {

///
/// Tests for Master-Master Connections
/// Acceptor and Requestor are different participants
///

template <typename T>
void TestSendAndReceivePrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", rank);
    {
      std::string msg("testOne");
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_TEST(msg == std::string("testTwo"));
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
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      std::string msg;
      com.receive(msg, 0);
      BOOST_TEST(msg == std::string("testOne"));
      msg = "testTwo";
      com.send(msg, 0);
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

template <typename T>
void TestSendAndReceiveVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", rank);
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
      std::vector<int> recv{1, 2, 3};
      com.receive(msg, 0);
      BOOST_TEST(msg == recv);
      com.send(msg, 0);
    }
    {
      std::vector<double> msg;
      com.receive(msg, 0);
      BOOST_TEST(msg == std::vector<double>({1.1, 2.2, 3.3}));
      com.send(msg, 0);
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
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
    com.closeConnection();
  }
}

template <typename T>
void TestSendReceiveFourProcesses(int rank)
{
  T   communication;
  int message = -1;

  switch (rank) {
  case 0: {
    communication.acceptConnection("A", "B", "", rank);

    communication.send(10, 2);
    communication.receive(message, 2);
    BOOST_TEST(message == 20);

    communication.send(20, 3);
    communication.receive(message, 3);
    BOOST_TEST(message == 40);

    communication.closeConnection();
    break;
  }
  case 1: {
    // does not participate
    break;
  }
  case 2: {
    communication.requestConnection("A", "B", "", rank, 2);

    communication.receive(message, 0);
    BOOST_TEST(message == 10);
    message *= 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  case 3: {
    communication.requestConnection("A", "B", "", rank, 2);

    communication.receive(message, 0);
    BOOST_TEST(message == 20);
    message *= 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  }
}

template <typename T>
void TestBroadcastPrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", rank);
    {
      double msg = 0.0;
      com.broadcast(msg);
    }
    {
      int msg = 1;
      com.broadcast(msg);
    }
    {
      bool msg = true;
      com.broadcast(msg);
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      double msg = 1.0;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == 0.0);
    }
    {
      int msg = 0;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == 1);
    }
    {
      bool msg = false;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == true);
    }
    com.closeConnection();
  }
}

template <typename T>
void TestBroadcastVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", rank);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(3.1415);
      com.broadcast(msg.data(), msg.size());
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(21);
      com.broadcast(msg.data(), msg.size());
    }
    {
      std::vector<int> msg{2, 3, 5, 8};
      com.broadcast(msg.data(), msg.size());
    }
    {
      std::vector<double> msg{1.2, 2.3, 3.5, 4.8};
      com.broadcast(msg.data(), msg.size());
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(3.1415)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(21)));
    }
    {
      std::vector<int> msg(4);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(msg == std::vector<int>({2, 3, 5, 8}));
    }
    {
      std::vector<double> msg(4);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(msg == std::vector<double>({1.2, 2.3, 3.5, 4.8}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReducePrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      int msg = 1;
      int rcv = 0;
      com.reduceSum(msg, rcv);
      BOOST_TEST(msg == 1);
      BOOST_TEST(rcv == 4);
    }
    {
      int msg = 1;
      int rcv = 0;
      com.allreduceSum(msg, rcv);
      BOOST_TEST(msg == 1);
      BOOST_TEST(rcv == 4);
    }
    {
      double msg = 3;
      double rcv = 0;
      com.allreduceSum(msg, rcv);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 3.1415);
    }

    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      int msg = 3;
      int rcv = 0;
      com.reduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 0);
    }
    {
      int msg = 3;
      int rcv = 0;
      com.allreduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 4);
    }
    {
      double msg = 0.1415;
      double rcv = 0;
      com.allreduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 0.1415);
      BOOST_TEST(rcv == 3.1415);
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReduceVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("process0", "process1", "", rank);
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg.data(), rcv.data(), msg.size());
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg.data(), rcv.data(), msg.size());
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg.data(), rcv.data(), msg.size(), 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{0, 0, 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg.data(), rcv.data(), msg.size(), 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  }
}

template <typename T>
void TestSendAndReceive(int rank)
{
  TestSendAndReceivePrimitiveTypes<T>(rank);
  TestSendAndReceiveVectors<T>(rank);
  TestBroadcastPrimitiveTypes<T>(rank);
  TestBroadcastVectors<T>(rank);
  TestReducePrimitiveTypes<T>(rank);
  TestReduceVectors<T>(rank);
}

} // namespace mastermaster

namespace masterslave {

///
/// Tests for Master-Slave Connections
/// Acceptor and Requestor are the same participant
///

template <typename T>
void TestSendAndReceivePrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", 0, 1);
    {
      std::string msg("testOne");
      com.send(msg, 1);
      com.receive(msg, 1);
      BOOST_TEST(msg == std::string("testTwo"));
    }
    {
      double msg = 0.0;
      com.send(msg, 1);
      com.receive(msg, 1);
      BOOST_TEST(msg == 1.0);
    }
    {
      int msg = 1;
      com.send(msg, 1);
      com.receive(msg, 1);
      BOOST_TEST(msg == 2);
    }
    {
      bool msg = true;
      com.send(msg, 1);
      com.receive(msg, 1);
      BOOST_TEST(msg == false);
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
    {
      std::string msg;
      com.receive(msg, 0);
      BOOST_TEST(msg == std::string("testOne"));
      msg = "testTwo";
      com.send(msg, 0);
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

template <typename T>
void TestSendAndReceiveVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.receive(msg.data(), msg.size(), 1);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(1)));
      msg = Eigen::Vector3d::Constant(2);
      com.send(msg.data(), msg.size(), 1);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.receive(msg.data(), msg.size(), 1);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(1)));
      msg = Eigen::Vector4i::Constant(0);
      com.send(msg.data(), msg.size(), 1);
    }
    {
      std::vector<int> msg;
      std::vector<int> recv{1, 2, 3};
      com.receive(msg, 1);
      BOOST_TEST(msg == recv);
      com.send(msg, 1);
    }
    {
      std::vector<double> msg;
      com.receive(msg, 1);
      BOOST_TEST(msg == std::vector<double>({1.1, 2.2, 3.3}));
      com.send(msg, 1);
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
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
    com.closeConnection();
  }
}

template <typename T>
void TestBroadcastPrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", rank, 1);
    {
      double msg = 0.0;
      com.broadcast(msg);
    }
    {
      int msg = 1;
      com.broadcast(msg);
    }
    {
      bool msg = true;
      com.broadcast(msg);
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
    {
      double msg = 1.0;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == 0.0);
    }
    {
      int msg = 0;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == 1);
    }
    {
      bool msg = false;
      com.broadcast(msg, 0);
      BOOST_TEST(msg == true);
    }
    com.closeConnection();
  }
}

template <typename T>
void TestBroadcastVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", rank, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(3.1415);
      com.broadcast(msg.data(), msg.size());
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(21);
      com.broadcast(msg.data(), msg.size());
    }
    {
      std::vector<int> msg{2, 3, 5, 8};
      com.broadcast(msg.data(), msg.size());
    }
    {
      std::vector<double> msg{1.2, 2.3, 3.5, 4.8};
      com.broadcast(msg.data(), msg.size());
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(3.1415)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(21)));
    }
    {
      std::vector<int> msg(4);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(msg == std::vector<int>({2, 3, 5, 8}));
    }
    {
      std::vector<double> msg(4);
      com.broadcast(msg.data(), msg.size(), 0);
      BOOST_CHECK(msg == std::vector<double>({1.2, 2.3, 3.5, 4.8}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReducePrimitiveTypes(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", 0, 1);
    {
      int msg = 1;
      int rcv = 0;
      com.reduceSum(msg, rcv);
      BOOST_TEST(msg == 1);
      BOOST_TEST(rcv == 4);
    }
    {
      int msg = 1;
      int rcv = 0;
      com.allreduceSum(msg, rcv);
      BOOST_TEST(msg == 1);
      BOOST_TEST(rcv == 4);
    }
    {
      double msg = 3;
      double rcv = 0;
      com.allreduceSum(msg, rcv);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 3.1415);
    }

    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
    {
      int msg = 3;
      int rcv = 0;
      com.reduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 0);
    }
    {
      int msg = 3;
      int rcv = 0;
      com.allreduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 3);
      BOOST_TEST(rcv == 4);
    }
    {
      double msg = 0.1415;
      double rcv = 0;
      com.allreduceSum(msg, rcv, 0);
      BOOST_TEST(msg == 0.1415);
      BOOST_TEST(rcv == 3.1415);
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReduceVectors(int rank)
{
  T com;

  if (rank == 0) {
    com.acceptConnection("Master", "Slave", "", rank, 1);
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg.data(), rcv.data(), msg.size());
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg.data(), rcv.data(), msg.size());
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  } else if (rank == 1) {
    com.requestConnection("Master", "Slave", "", 0, 1);
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg.data(), rcv.data(), msg.size(), 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{0, 0, 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg.data(), rcv.data(), msg.size(), 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(), msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(), rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  }
}

template <typename T>
void TestSendAndReceive(int rank)
{
  TestSendAndReceivePrimitiveTypes<T>(rank);
  TestSendAndReceiveVectors<T>(rank);
  TestBroadcastPrimitiveTypes<T>(rank);
  TestBroadcastVectors<T>(rank);
  TestReducePrimitiveTypes<T>(rank);
  TestReduceVectors<T>(rank);
}

} // namespace masterslave

namespace serverclient {

///
/// Tests for Server-Client Connections
/// Server awaits n connections
/// n Clients connect to the server
/// The Server and the set of Clients are on different participants

/// Tests connecting two processes using acceptConnectionAsServer and requestConnectionAsClient
template <typename T>
void TestSendReceiveTwoProcessesServerClient(int rank)
{
  T   communication;
  int message = 1;

  switch (rank) {
  case 0: {
    communication.acceptConnectionAsServer("A", "B", "", rank, 1);
    communication.send(message, 1);
    communication.receive(message, 1);
    BOOST_TEST(message == 2);
    communication.closeConnection();
    break;
  }
  case 1: {
    communication.requestConnectionAsClient("A", "B", "", {0}, rank);
    communication.receive(message, 0);
    BOOST_TEST(message == 1);
    message = 2;
    communication.send(message, 0);
    communication.closeConnection();
    break;
  }
  }
}

template <typename T>
void TestSendReceiveFourProcessesServerClient(int rank)
{
  T   communication;
  int message = -1;

  switch (rank) {
  case 0: {
    communication.acceptConnectionAsServer("A", "B", "", rank, 2);

    communication.send(10, 2);
    communication.receive(message, 2);
    BOOST_TEST(message == 20);

    communication.send(20, 3);
    communication.receive(message, 3);
    BOOST_TEST(message == 40);

    communication.closeConnection();
    break;
  }
  case 1: {
    // does not accept a connection
    break;
  }
  case 2: {
    communication.requestConnectionAsClient("A", "B", "", {0}, rank);

    communication.receive(message, 0);
    BOOST_TEST(message == 10);
    message *= 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  case 3: {
    communication.requestConnectionAsClient("A", "B", "", {0}, rank);

    communication.receive(message, 0);
    BOOST_TEST(message == 20);
    message *= 2;
    communication.send(message, 0);

    communication.closeConnection();
    break;
  }
  }
}

template <typename T>
void TestSendReceiveFourProcessesServerClientV2(int rank)
{
  T   communication;
  int message = -1;

  switch (rank) {
  case 0: {
    communication.acceptConnectionAsServer("A", "B", "", rank, 2);

    communication.send(10, 2);
    communication.receive(message, 2);
    BOOST_TEST(message == 20);

    communication.send(100, 3);
    communication.receive(message, 3);
    BOOST_TEST(message == 200);

    communication.closeConnection();
    break;
  }
  case 1: {
    communication.acceptConnectionAsServer("A", "B", "", rank, 2);

    communication.send(20, 2);
    communication.receive(message, 2);
    BOOST_TEST(message == 40);

    communication.send(200, 3);
    communication.receive(message, 3);
    BOOST_TEST(message == 400);

    communication.closeConnection();
    break;
  }
  case 2: {
    communication.requestConnectionAsClient("A", "B", "", {0, 1}, rank);

    communication.receive(message, 0);
    BOOST_TEST(message == 10);
    communication.send(20, 0);

    communication.receive(message, 1);
    BOOST_TEST(message == 20);
    communication.send(40, 1);

    communication.closeConnection();
    break;
  }
  case 3: {
    communication.requestConnectionAsClient("A", "B", "", {0, 1}, rank);

    communication.receive(message, 0);
    BOOST_TEST(message == 100);
    communication.send(200, 0);

    communication.receive(message, 1);
    BOOST_TEST(message == 200);
    communication.send(400, 1);

    communication.closeConnection();
    break;
  }
  }
}

} // namespace serverclient

} // namespace com
} // namespace testing
} // namespace precice

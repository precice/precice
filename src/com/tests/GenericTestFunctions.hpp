#pragma once

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <string>
#include <vector>

#include "com/Communication.hpp"
#include "testing/Testing.hpp"

/// Generic test function that is called from the tests for
/// MPIPortsCommunication, MPIDirectCommunication and SocketCommunication

namespace precice {
namespace testing {
namespace com {

namespace primaryprimary {

///
/// Tests for primary connections
/// Acceptor and Requestor are different participants
///

template <typename T>
void TestSendAndReceivePrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
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
  } else {
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
void TestSendAndReceiveEigen(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(1)));
      msg = Eigen::Vector3d::Constant(2);
      com.send(msg, 0);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(1)));
      msg = Eigen::Vector4i::Constant(0);
      com.send(msg, 0);
    }
    com.closeConnection();
  } else {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(1);
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(2)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(1);
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Zero()));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestSendAndReceiveRanges(TestContext const &context)
{
  using precice::com::AsVectorTag;
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      std::vector<int> recv{1, 2, 3};
      std::vector<int> msg = com.receiveRange(0, AsVectorTag<int>{});
      BOOST_TEST(msg == recv);
      com.sendRange(msg, 0);
    }
    {
      std::vector<double> msg = com.receiveRange(0, AsVectorTag<double>{});
      BOOST_TEST(msg == std::vector<double>({1.1, 2.2, 3.3}));
      com.sendRange(msg, 0);
    }
    com.closeConnection();
  } else {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      std::vector<int> msg{1, 2, 3};
      com.sendRange(msg, 0);
      msg = com.receiveRange(0, AsVectorTag<int>{});
      BOOST_CHECK(msg == std::vector<int>({1, 2, 3}));
    }
    {
      std::vector<double> msg{1.1, 2.2, 3.3};
      com.sendRange(msg, 0);
      msg = com.receiveRange(0, AsVectorTag<double>{});
      BOOST_CHECK(msg == std::vector<double>({1.1, 2.2, 3.3}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestSendReceiveFourProcesses(TestContext const &context)
{
  T   communication;
  int message = -1;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      communication.acceptConnection("A", "B", "", 0);

      communication.send(10, 0);
      communication.receive(message, 0);
      BOOST_TEST(message == 20);

      communication.send(20, 1);
      communication.receive(message, 1);
      BOOST_TEST(message == 40);

      communication.closeConnection();
    }
  } else {
    if (context.isPrimary()) {
      communication.requestConnection("A", "B", "", 0, 2);

      communication.receive(message, 0);
      BOOST_TEST(message == 10);
      message *= 2;
      communication.send(message, 0);

      communication.closeConnection();
    } else {
      communication.requestConnection("A", "B", "", 1, 2);

      communication.receive(message, 0);
      BOOST_TEST(message == 20);
      message *= 2;
      communication.send(message, 0);

      communication.closeConnection();
    }
  }
}

template <typename T>
void TestBroadcastPrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
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
  } else {
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
void TestBroadcastEigen(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(3.1415);
      com.broadcast(msg);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(21);
      com.broadcast(msg);
    }
    com.closeConnection();
  } else {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.broadcast(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(3.1415)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.broadcast(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(21)));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestBroadcastVectors(TestContext const &context)
{
  T com;
  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      std::vector<int> msg{2, 3, 5, 8};
      com.broadcast(msg);
    }
    {
      std::vector<double> msg{1.2, 2.3, 3.5, 4.8};
      com.broadcast(msg);
    }
    com.closeConnection();
  } else {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      std::vector<int> msg(4);
      com.broadcast(msg, 0);
      BOOST_CHECK(msg == std::vector<int>({2, 3, 5, 8}));
    }
    {
      std::vector<double> msg(4);
      com.broadcast(msg, 0);
      BOOST_CHECK(msg == std::vector<double>({1.2, 2.3, 3.5, 4.8}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReducePrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
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
  } else {
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
void TestReduceVectors(TestContext const &context)
{
  T com;

  if (context.isNamed("A")) {
    com.acceptConnection("process0", "process1", "", 0);
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg, rcv);
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg, rcv);
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  } else {
    com.requestConnection("process0", "process1", "", 0, 1);
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg, rcv, 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{0, 0, 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg, rcv, 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  }
}

} // namespace primaryprimary

namespace intracomm {

///
/// Tests for intra-participant communication Connections
/// Acceptor and Requestor are the same participant
///

template <typename T>
void TestSendAndReceivePrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
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
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
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
void TestSendAndReceiveEigen(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.receive(msg, 1);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(1)));
      msg = Eigen::Vector3d::Constant(2);
      com.send(msg, 1);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.receive(msg, 1);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(1)));
      msg = Eigen::Vector4i::Constant(0);
      com.send(msg, 1);
    }
    com.closeConnection();
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(1);
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(2)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(1);
      com.send(msg, 0);
      com.receive(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Zero()));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestSendAndReceiveRanges(TestContext const &context)
{
  T com;
  using precice::com::AsVectorTag;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
    {
      std::vector<int> recv{1, 2, 3};
      std::vector<int> msg = com.receiveRange(1, AsVectorTag<int>{});
      BOOST_TEST(msg == recv);
      com.sendRange(msg, 1);
    }
    {
      std::vector<double> msg = com.receiveRange(1, AsVectorTag<double>{});
      BOOST_TEST(msg == std::vector<double>({1.1, 2.2, 3.3}));
      com.sendRange(msg, 1);
    }
    com.closeConnection();
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
    {
      std::vector<int> msg{1, 2, 3};
      com.sendRange(msg, 0);
      msg = com.receiveRange(0, AsVectorTag<int>{});
      BOOST_CHECK(msg == std::vector<int>({1, 2, 3}));
    }
    {
      std::vector<double> msg{1.1, 2.2, 3.3};
      com.sendRange(msg, 0);
      msg = com.receiveRange(0, AsVectorTag<double>{});
      BOOST_CHECK(msg == std::vector<double>({1.1, 2.2, 3.3}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestBroadcastPrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
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
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
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
void TestBroadcastVectors(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(3.1415);
      com.broadcast(msg);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(21);
      com.broadcast(msg);
    }
    {
      std::vector<int> msg{2, 3, 5, 8};
      com.broadcast(msg);
    }
    {
      std::vector<double> msg{1.2, 2.3, 3.5, 4.8};
      com.broadcast(msg);
    }
    com.closeConnection();
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.broadcast(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector3d::Constant(3.1415)));
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.broadcast(msg, 0);
      BOOST_CHECK(testing::equals(msg, Eigen::Vector4i::Constant(21)));
    }
    {
      std::vector<int> msg(4);
      com.broadcast(msg, 0);
      BOOST_CHECK(msg == std::vector<int>({2, 3, 5, 8}));
    }
    {
      std::vector<double> msg(4);
      com.broadcast(msg, 0);
      BOOST_CHECK(msg == std::vector<double>({1.2, 2.3, 3.5, 4.8}));
    }
    com.closeConnection();
  }
}

template <typename T>
void TestReducePrimitiveTypes(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
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
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
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
void TestReduceVectors(TestContext const &context)
{
  T com;

  if (context.isPrimary()) {
    com.acceptConnection("Primary", "Secondary", "", 0, 1);
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg, rcv);
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{0.1, 0.2, 0.3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg, rcv);
      std::vector<double> msg_expected{0.1, 0.2, 0.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  } else {
    com.requestConnection("Primary", "Secondary", "", 0, 1);
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.reduceSum(msg, rcv, 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{0, 0, 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    {
      std::vector<double> msg{1, 2, 3};
      std::vector<double> rcv{0, 0, 0};
      com.allreduceSum(msg, rcv, 0);
      std::vector<double> msg_expected{1, 2, 3};
      BOOST_CHECK_EQUAL_COLLECTIONS(msg.begin(), msg.end(),
                                    msg_expected.begin(), msg_expected.end());
      std::vector<double> rcv_expected{1.1, 2.2, 3.3};
      BOOST_CHECK_EQUAL_COLLECTIONS(rcv.begin(), rcv.end(),
                                    rcv_expected.begin(), rcv_expected.end());
    }
    com.closeConnection();
  }
}

} // namespace intracomm

namespace serverclient {

///
/// Tests for Server-Client Connections
/// Server awaits n connections
/// n Clients connect to the server
/// The Server and the set of Clients are on different participants

/// Tests connecting two processes using acceptConnectionAsServer and
/// requestConnectionAsClient
template <typename T>
void TestSendReceiveTwoProcessesServerClient(TestContext const &context)
{
  T   communication;
  int message = 1;
  BOOST_REQUIRE(context.hasSize(1));

  if (context.isNamed("A")) {
    communication.acceptConnectionAsServer("A", "B", "", 0, 1);
    communication.send(message, 0);
    communication.receive(message, 0);
    BOOST_TEST(message == 2);
    communication.closeConnection();
  } else {
    BOOST_REQUIRE(context.isNamed("B"));
    communication.requestConnectionAsClient("A", "B", "", {0}, 0);
    communication.receive(message, 0);
    BOOST_TEST(message == 1);
    message = 2;
    communication.send(message, 0);
    communication.closeConnection();
  }
}

template <typename T>
void TestSendReceiveFourProcessesServerClient(TestContext const &context)
{
  T   communication;
  int message = -1;
  BOOST_REQUIRE(context.hasSize(2));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      communication.acceptConnectionAsServer("A", "B", "", context.rank, 2);

      communication.send(10, 0);
      communication.receive(message, 0);
      BOOST_TEST(message == 20);

      communication.send(20, 1);
      communication.receive(message, 1);
      BOOST_TEST(message == 40);

      communication.closeConnection();
    } else {
      communication.acceptConnectionAsServer("A", "B", "", context.rank, 0);
      communication.closeConnection();
    }
  } else {
    BOOST_REQUIRE(context.isNamed("B"));
    if (context.isPrimary()) {
      communication.requestConnectionAsClient("A", "B", "", {0}, context.rank);

      communication.receive(message, 0);
      BOOST_TEST(message == 10);
      message *= 2;
      communication.send(message, 0);

      communication.closeConnection();
    } else {
      communication.requestConnectionAsClient("A", "B", "", {0}, context.rank);

      communication.receive(message, 0);
      BOOST_TEST(message == 20);
      message *= 2;
      communication.send(message, 0);

      communication.closeConnection();
    }
  }
}

template <typename T>
void TestSendReceiveFourProcessesServerClientV2(TestContext const &context)
{
  T   communication;
  int message = -1;

  BOOST_REQUIRE(context.hasSize(2));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      communication.acceptConnectionAsServer("A", "B", "", 0, 2);

      communication.send(10, 0);
      communication.receive(message, 0);
      BOOST_TEST(message == 20);

      communication.send(100, 1);
      communication.receive(message, 1);
      BOOST_TEST(message == 200);

      communication.closeConnection();
    } else {
      communication.acceptConnectionAsServer("A", "B", "", 1, 2);

      communication.send(20, 0);
      communication.receive(message, 0);
      BOOST_TEST(message == 40);

      communication.send(200, 1);
      communication.receive(message, 1);
      BOOST_TEST(message == 400);

      communication.closeConnection();
    }

  } else {
    BOOST_REQUIRE(context.isNamed("B"));

    if (context.isPrimary()) {
      communication.requestConnectionAsClient("A", "B", "", {0, 1}, 0);

      communication.receive(message, 0);
      BOOST_TEST(message == 10);
      communication.send(20, 0);

      communication.receive(message, 1);
      BOOST_TEST(message == 20);
      communication.send(40, 1);

      communication.closeConnection();
    } else {
      communication.requestConnectionAsClient("A", "B", "", {0, 1}, 1);

      communication.receive(message, 0);
      BOOST_TEST(message == 100);
      communication.send(200, 0);

      communication.receive(message, 1);
      BOOST_TEST(message == 200);
      communication.send(400, 1);

      communication.closeConnection();
    }
  }
}

} // namespace serverclient

} // namespace com
} // namespace testing
} // namespace precice

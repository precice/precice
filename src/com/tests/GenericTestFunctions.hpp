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

/// This tests still uses the old rank enumeration
template <typename T>
void TestSendAndReceivePrimitiveTypes()
{
  T         com;
  const int rank = utils::Parallel::getProcessRank();

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

/// This tests still uses the old rank enumeration
template <typename T>
void TestSendAndReceiveVectors()
{
  T         com;
  const int rank = utils::Parallel::getProcessRank();

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

/// Tests connecting four processes using acceptConnection and requestConnection
template <typename T>
void TestSendReceiveFourProcesses()
{
  T         communication;
  const int rank    = utils::Parallel::getProcessRank();
  int       message = -1;

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
void TestSendAndReceive()
{
  TestSendAndReceivePrimitiveTypes<T>();
  TestSendAndReceiveVectors<T>();
}

/// Tests connecting two processes using acceptConnectionAsServer and requestConnectionAsClient
template <typename T>
void TestSendReceiveTwoProcessesServerClient()
{
  T         communication;
  const int rank    = utils::Parallel::getProcessRank();
  int       message = 1;

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
void TestSendReceiveFourProcessesServerClient()
{
  T         communication;
  const int rank    = utils::Parallel::getProcessRank();
  int       message = -1;

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
void TestSendReceiveFourProcessesServerClientV2()
{
  T   communication;
  int rank    = utils::Parallel::getProcessRank();
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

#include <memory>
#include <string>
#include <vector>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_SUITE(Parallel)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(RestrictCommTest)
{
  PRECICE_TEST(4_ranks);
  using Par = utils::Parallel;

  auto baseComm = Par::current();
  Par::restrictCommunicator(3);
  auto restComm = Par::current();
  BOOST_TEST(restComm);
  BOOST_TEST(restComm->parent == baseComm);

  if (context.rank == 3) {
    BOOST_TEST(restComm->isNull());
  } else {
    BOOST_TEST(!restComm->isNull());
    BOOST_TEST(restComm->size() == 3);
    BOOST_TEST(restComm->rank() <= 2);
  }
}

BOOST_AUTO_TEST_CASE(SplitCommTest)
{
  PRECICE_TEST(3_ranks);
  using Par = utils::Parallel;

  std::string name{"Group"};
  name.append(context.rank == 2 ? "Two" : "One");
  auto baseComm = Par::current();
  Par::splitCommunicator(name);
  auto splitComm = Par::current();
  BOOST_TEST(splitComm->parent == baseComm);

  BOOST_TEST(!splitComm->isNull());

  const auto &groups = splitComm->groups;
  BOOST_TEST(groups.size() == 2);
  BOOST_TEST(groups.at(0).id == 0);
  BOOST_TEST(groups.at(1).id == 1);
  BOOST_TEST(groups.at(0).name == std::string("GroupOne"));
  BOOST_TEST(groups.at(1).name == std::string("GroupTwo"));
  BOOST_TEST(groups.at(0).leaderRank == 0);
  BOOST_TEST(groups.at(1).leaderRank == 2);
  BOOST_TEST(groups.at(0).size == 2);
  BOOST_TEST(groups.at(1).size == 1);
}

BOOST_AUTO_TEST_CASE(Master1SlaveTest)
{
  PRECICE_TEST(""_on(2_ranks).setupMasterSlaves());

  BOOST_TEST(context.hasSize(2));
  auto &com = precice::utils::MasterSlave::_communication;
  BOOST_TEST((com != nullptr));

  if (context.isMaster()) {
    int first = 1001;
    com->send(first, 1);

    int answer = -1;
    com->receive(answer, 1);
    BOOST_TEST(answer == 1111);
  } else {
    int received = -1;
    com->receive(received, 0);
    BOOST_TEST(received == 1001);

    received += 110;
    BOOST_TEST(received == 1111);
    com->send(received, 0);
  }
}

BOOST_AUTO_TEST_CASE(Master2SlaveTest)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  BOOST_TEST(context.hasSize(3));
  auto &com = precice::utils::MasterSlave::_communication;
  BOOST_TEST((com != nullptr));

  if (context.isMaster()) {
    int sum = 1001;

    int message = -1;
    com->receive(message, 1);
    sum += message;
    BOOST_TEST(sum == 1011);

    message = -1;
    com->receive(message, 2);
    sum += message;
    BOOST_TEST(sum == 1111);

    com->send(sum, 1);
    com->send(sum, 2);
  } else {
    int tosend = context.isRank(1) ? 10 : 100;
    com->send(tosend, 0);

    int received = -1;
    com->receive(received, 0);
    BOOST_TEST(received == 1111);
  }
}

BOOST_AUTO_TEST_CASE(OffsetMaster1SlaveTest)
{
  PRECICE_TEST("Offset"_on(1_rank), "Test"_on(2_ranks).setupMasterSlaves());

  if (context.isNamed("Offset"))
    return;

  BOOST_TEST(context.hasSize(2));
  auto &com = precice::utils::MasterSlave::_communication;
  BOOST_TEST((com != nullptr));

  if (context.isMaster()) {
    int first = 1001;

    com->send(first, 1);

    int answer = -1;
    com->receive(answer, 1);
    BOOST_TEST(answer == 1111);
  } else {
    int received = -1;
    com->receive(received, 0);
    BOOST_TEST(received == 1001);

    received += 110;
    com->send(received, 0);
  }
}

BOOST_AUTO_TEST_CASE(OffsetMaster2SlaveTest)
{
  PRECICE_TEST("Offset"_on(1_rank), "Test"_on(3_ranks).setupMasterSlaves());

  if (context.isNamed("Offset"))
    return;

  BOOST_TEST(context.hasSize(3));
  auto &com = precice::utils::MasterSlave::_communication;
  BOOST_TEST((com != nullptr));

  if (context.isMaster()) {
    int sum = 1001;

    int message = -1;
    com->receive(message, 1);
    sum += message;
    BOOST_TEST(sum == 1011);

    message = -1;
    com->receive(message, 2);
    sum += message;
    BOOST_TEST(sum == 1111);

    com->send(sum, 1);
    com->send(sum, 2);
  } else {
    int tosend = context.isRank(1) ? 10 : 100;
    com->send(tosend, 0);

    int received = -1;
    com->receive(received, 0);
    BOOST_TEST(received == 1111);
  }
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

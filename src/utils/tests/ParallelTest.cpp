#include <memory>
#include <string>
#include <vector>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_SUITE(Parallel)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(Primary1SecondaryTest)
{
  PRECICE_TEST(""_on(2_ranks).setupIntraComm());

  BOOST_TEST(context.hasSize(2));
  auto &com = precice::utils::IntraComm::getCommunication();
  BOOST_TEST((com != nullptr));

  if (context.isPrimary()) {
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

BOOST_AUTO_TEST_CASE(Primary2SecondaryTest)
{
  PRECICE_TEST(""_on(3_ranks).setupIntraComm());

  BOOST_TEST(context.hasSize(3));
  auto &com = precice::utils::IntraComm::getCommunication();
  BOOST_TEST((com != nullptr));

  if (context.isPrimary()) {
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

BOOST_AUTO_TEST_CASE(OffsetPrimary1SecondaryTest)
{
  PRECICE_TEST("Offset"_on(1_rank), "Test"_on(2_ranks).setupIntraComm());

  if (context.isNamed("Offset"))
    return;

  BOOST_TEST(context.hasSize(2));
  auto &com = precice::utils::IntraComm::getCommunication();
  BOOST_TEST((com != nullptr));

  if (context.isPrimary()) {
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

BOOST_AUTO_TEST_CASE(OffsetPrimary2SecondaryTest)
{
  PRECICE_TEST("Offset"_on(1_rank), "Test"_on(3_ranks).setupIntraComm());

  if (context.isNamed("Offset"))
    return;

  BOOST_TEST(context.hasSize(3));
  auto &com = precice::utils::IntraComm::getCommunication();
  BOOST_TEST((com != nullptr));

  if (context.isPrimary()) {
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

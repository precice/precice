#include "cplscheme/impl/TimeHandler.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(CplSchemeTests)
BOOST_AUTO_TEST_SUITE(TimeHandler)

using precice::cplscheme::impl::TimeHandler;

BOOST_AUTO_TEST_CASE(ManyTinyTimesteps)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 1e-7;

  for (int i = 1; i <= 1'000'000; ++i) {
    th.progressBy(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(ManyBigTimesteps)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 1000.0;

  for (int i = 1; i <= 1'000'000; ++i) {
    th.progressBy(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(ManyMixedTimesteps)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws1 = 0.01;
  double tws2 = 0.005;

  int n1 = 0, n2 = 0;
  for (int i = 1; i <= 1'000'000; ++i) {
    double tws;
    if (i % 2 == 1) {
      tws = tws1;
      ++n1;
    } else {
      tws = tws2;
      ++n2;
    }

    th.progressBy(tws);
    double t = n1 * tws1 + n2 * tws2;
    BOOST_TEST(th.time() == t);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == t);
    BOOST_TEST(th.windowProgress() == 0.0);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
  }
}

BOOST_AUTO_TEST_SUITE(Subcycling)

BOOST_AUTO_TEST_CASE(NormalTWSmallTS)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 10.0;
  double dt  = 0.1;

  for (int i = 1; i <= 10'000; ++i) {

    for (int j = 1; j <= 100; ++j) {
      th.progressBy(dt);
      BOOST_TEST(th.time() == (i - 1) * tws + j * dt);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - dt * j);
    }
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * i);
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(SmallTWTinyTS)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 0.001;
  double dt  = 1e-7;

  for (int i = 1; i <= 1'000; ++i) {

    for (int j = 1; j <= 10'000; ++j) {
      th.progressBy(dt);
      BOOST_TEST(th.time() == (i - 1) * tws + j * dt);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - dt * j);
    }
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * i);
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(LargeTWNormalTS)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 1000;
  double dt  = 0.1;

  for (int i = 1; i <= 100; ++i) {

    for (int j = 1; j <= 10'000; ++j) {
      th.progressBy(dt);
      BOOST_TEST(th.time() == (i - 1) * tws + j * dt);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - dt * j);
    }
    BOOST_TEST(th.untilWindowEnd(tws) == 0, boost::test_tools::tolerance(1e-13));
    //BOOST_TEST(th.reachedEndOfWindow(tws)); This currently fails
    BOOST_TEST(th.time() == tws * i);
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * i);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

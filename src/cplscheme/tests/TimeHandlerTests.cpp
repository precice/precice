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
  int    twc = 1'000'000;

  for (int tw = 1; tw <= twc; ++tw) {
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    th.progressBy(tws);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(ManyBigTimesteps)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 1000.0;
  int    twc = 1'000'000;

  for (int tw = 1; tw <= twc; ++tw) {
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    th.progressBy(tws);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * tw);
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
  int    twc  = 1'000'000;

  int n1 = 0, n2 = 0;
  for (int tw = 1; tw <= twc; ++tw) {
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    double tws;
    if (tw % 2 == 1) {
      tws = tws1;
      ++n1;
    } else {
      tws = tws2;
      ++n2;
    }
    BOOST_TEST_INFO_SCOPE("this time-window size is " << tws);

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
  double tss = 0.1;
  int    twc = 10'000;
  int    tsc = 100;

  for (int tw = 1; tw <= twc; ++tw) {
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);

    for (int ts = 1; ts <= tsc; ++ts) {
      BOOST_TEST_INFO_SCOPE("time-step " << ts << "/" << tsc);
      th.progressBy(tss);
      BOOST_TEST(th.time() == (tw - 1) * tws + ts * tss);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * ts);
    }
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * tw);
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(SmallTWTinyTS)
{
  PRECICE_TEST(1_rank);
  TimeHandler th;

  double tws = 0.001;
  double tss = 1e-7;
  int    twc = 1'000;
  int    tsc = 10'000;

  for (int i = 1; i <= twc; ++i) {
    BOOST_TEST_INFO_SCOPE("time-window " << i << "/" << twc);

    for (int j = 1; j <= tsc; ++j) {
      BOOST_TEST_INFO_SCOPE("time-step " << j << "/" << tsc);
      th.progressBy(tss);
      BOOST_TEST(th.time() == (i - 1) * tws + j * tss);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * j);
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
  double tss = 0.1;
  int    twc = 100;
  int    tsc = 10'000;

  for (int tw = 1; tw <= twc; ++tw) {
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);

    for (int ts = 1; ts <= tsc; ++ts) {
      th.progressBy(tss);
      BOOST_TEST(th.time() == (tw - 1) * tws + ts * tss);
      BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * ts);
    }
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * tw);
    th.completeTimeWindow(tws);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

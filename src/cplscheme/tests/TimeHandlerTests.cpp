#include "cplscheme/impl/TimeHandler.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(CplSchemeTests)
BOOST_AUTO_TEST_SUITE(TimeHandler)

using precice::cplscheme::impl::TimeHandler;

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ManyTinyTimesteps)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws        = 1e-7;
  int    twc        = 1'000'000;
  int    checkEvery = 1000;

  for (int tw = 1; tw <= twc; ++tw) {
    th.progressBy(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.time() == tws * tw);
      BOOST_TEST(th.reachedEndOfWindow(tws));
    }
    th.completeTimeWindow(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.time() == tws * tw);
      BOOST_TEST(!th.reachedEndOfWindow(tws));
      BOOST_TEST(th.windowProgress() == 0.0);
    }
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ManyBigTimesteps)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws        = 1000.0;
  int    twc        = 1'000'000;
  int    checkEvery = 1000;

  for (int tw = 1; tw <= twc; ++tw) {
    th.progressBy(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.time() == tws * tw);
      BOOST_TEST(th.reachedEndOfWindow(tws));
    }
    th.completeTimeWindow(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.time() == tws * tw);
      BOOST_TEST(!th.reachedEndOfWindow(tws));
      BOOST_TEST(th.windowProgress() == 0.0);
    }
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ManyMixedTimesteps)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws1       = 0.01;
  double tws2       = 0.005;
  int    twc        = 1'000'000;
  int    checkEvery = 1000;

  int n1 = 0, n2 = 0;
  for (int tw = 1; tw <= twc; ++tw) {
    double tws;
    if (tw % 2 == 1) {
      tws = tws1;
      ++n1;
    } else {
      tws = tws2;
      ++n2;
    }
    th.progressBy(tws);
    double t = n1 * tws1 + n2 * tws2;
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST_INFO_SCOPE("this time-window size is " << tws);
      BOOST_TEST(th.time() == t);
      BOOST_TEST(th.reachedEndOfWindow(tws));
    }
    th.completeTimeWindow(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST_INFO_SCOPE("this time-window size is " << tws);
      BOOST_TEST(th.time() == t);
      BOOST_TEST(th.windowProgress() == 0.0);
      BOOST_TEST(!th.reachedEndOfWindow(tws));
    }
  }
}

BOOST_AUTO_TEST_SUITE(Subcycling)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NormalTWSmallTS)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws        = 10.0;
  double tss        = 0.1;
  int    twc        = 10'000;
  int    tsc        = 100;
  int    checkEvery = 1000;

  for (int tw = 1; tw <= twc; ++tw) {

    for (int ts = 1; ts <= tsc; ++ts) {
      th.progressBy(tss);
      if (tw % checkEvery == 0) {
        BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
        BOOST_TEST_INFO_SCOPE("time-step " << ts << "/" << tsc);
        BOOST_TEST(th.time() == (tw - 1) * tws + ts * tss);
        BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * ts);
      }
    }
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.reachedEndOfWindow(tws));
      BOOST_TEST(th.time() == tws * tw);
    }
    th.completeTimeWindow(tws);
    if (tw % checkEvery == 0) {
      BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
      BOOST_TEST(th.time() == tws * tw);
      BOOST_TEST(!th.reachedEndOfWindow(tws));
      BOOST_TEST(th.windowProgress() == 0.0);
    }
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SmallTWTinyTS)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws        = 0.001;
  double tss        = 1e-7;
  int    twc        = 100;
  int    tsc        = 10'000;
  int    checkEvery = 1000;

  for (int tw = 1; tw <= twc; ++tw) {
    for (int ts = 1; ts <= tsc; ++ts) {
      th.progressBy(tss);
      if (ts % checkEvery == 0) {
        BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
        BOOST_TEST_INFO_SCOPE("time-step " << ts << "/" << tsc);
        BOOST_TEST(th.time() == (tw - 1) * tws + ts * tss);
        BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * ts);
      }
    }
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * tw);
    th.completeTimeWindow(tws);
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(LargeTWNormalTS)
{
  PRECICE_TEST();
  TimeHandler th;

  double tws        = 1000;
  double tss        = 0.1;
  int    twc        = 100;
  int    tsc        = 10'000;
  int    checkEvery = 1000;

  for (int tw = 1; tw <= twc; ++tw) {
    for (int ts = 1; ts <= tsc; ++ts) {
      th.progressBy(tss);
      if (ts % checkEvery == 0) {
        BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
        BOOST_TEST_INFO_SCOPE("time-step " << ts << "/" << tsc);
        BOOST_TEST(th.time() == (tw - 1) * tws + ts * tss);
        BOOST_TEST(th.untilWindowEnd(tws) == tws - tss * ts);
      }
    }
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    BOOST_TEST(th.reachedEndOfWindow(tws));
    BOOST_TEST(th.time() == tws * tw);
    th.completeTimeWindow(tws);
    BOOST_TEST_INFO_SCOPE("time-window " << tw << "/" << twc);
    BOOST_TEST(th.time() == tws * tw);
    BOOST_TEST(!th.reachedEndOfWindow(tws));
    BOOST_TEST(th.windowProgress() == 0.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MaxTime)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(StepToEnd)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  th.progressBy(1.0);
  BOOST_TEST(th.time() == 1.0);
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);

  th.completeTimeWindow(1.0);
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(StepTruncated)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tws = 2.0;

  th.progressBy(1.0);

  BOOST_TEST(th.time() == 1.0);
  BOOST_TEST(th.reachedEndOfWindow(tws));
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.untilWindowEnd(tws) == 0.0);

  th.completeTimeWindow(tws);
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(tws));
  BOOST_TEST(th.untilWindowEnd(tws) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(StepCloseToEnd)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  th.progressBy(0.999999999999999999999999);
  BOOST_TEST(th.time() == 1.0);
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);

  th.completeTimeWindow(1.0);
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TinyTimestepsToEnd)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tss = 1e-6;
  double tsc = 1'000'000;

  for (int ts = 1; ts <= tsc; ++ts) {
    th.progressBy(tss);
  }

  BOOST_TEST(th.time() == 1.0);
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);

  th.completeTimeWindow(1.0);
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ManyTimeWindowsToEnd)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tws = 1e-6;
  double twc = 1'000'000;

  for (int tw = 1; tw <= twc; ++tw) {
    th.progressBy(tws);
    th.completeTimeWindow(tws);
  }

  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TimeWindowsWithSubstepsToEnd)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tws = 1e-3;
  double twc = 1'000;
  double tss = 1e-6;
  double tsc = 1'000;

  for (int tw = 1; tw <= twc; ++tw) {
    for (int ts = 1; ts <= tsc; ++ts) {
      th.progressBy(tss);
    }
    BOOST_TEST(th.reachedEndOfWindow(tws));
    th.completeTimeWindow(tws);
  }
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TinyTimestepsTruncated)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tss = 1e-6;
  double tsc = 1'000'000;

  for (int ts = 1; ts <= tsc; ++ts) {
    th.progressBy(tss);
  }

  BOOST_TEST(th.time() == 1.0);
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);

  th.completeTimeWindow(1.0);
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ManyTimeWindowsTrucated)
{
  PRECICE_TEST();
  TimeHandler th(1.0);

  double tws = 1e-6;
  double twc = 1'000'000;
  double twl = 0.001;

  for (int tw = 1; tw <= twc; ++tw) {
    th.progressBy(tws);
    if (tw == twc) {
      // extend last time window past the end, so it is truncated my maxTime
      th.completeTimeWindow(twl);
    } else {
      th.completeTimeWindow(tws);
    }
  }
  BOOST_TEST(th.windowProgress() == 0.0);
  BOOST_TEST(th.reachedEnd());
  BOOST_TEST(th.reachedEndOfWindow(1.0));
  BOOST_TEST(th.untilWindowEnd(1.0) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

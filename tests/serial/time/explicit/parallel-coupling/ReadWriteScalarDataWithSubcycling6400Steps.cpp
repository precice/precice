#ifndef PRECICE_NO_MPI

#include "math/differences.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>
#include "helpers.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * This is a smoke test to reproduce the scenario explained in https://github.com/precice/precice/issues/1866
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcycling6400Steps)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  bool useAdvancedDtStrategy = true;

  if (context.isNamed("SolverOne")) {
    const int              nSubsteps = 6400;
    const std::vector<int> expectedSubsteps{nSubsteps,
                                            (useAdvancedDtStrategy) ? nSubsteps + 1 : nSubsteps + 1, // @todo performing 6401 instead of 6400 steps in second window due to round-off errors inside preCICE, see https://github.com/precice/precice/issues/1866. At least with useAdvancedDtStrategy = true, a user should be able to work-around this problem.
                                            nSubsteps, nSubsteps, nSubsteps};
    subcyclingWithNSteps(context, nSubsteps, expectedSubsteps, useAdvancedDtStrategy);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    const int nSubsteps = 1;
    subcyclingWithNSteps(context, nSubsteps, std::vector<int>{nSubsteps, nSubsteps, nSubsteps, nSubsteps, nSubsteps}, useAdvancedDtStrategy);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

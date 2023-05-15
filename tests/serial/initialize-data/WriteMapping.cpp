#ifndef PRECICE_NO_MPI

#include <precice/Participant.hpp>
#include <vector>
#include "helpers.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(InitializeData)

/**
 * @brief The second solver initializes the data of the first. Use write mapping for data.
 */
BOOST_AUTO_TEST_CASE(WriteMapping)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  testDataInitialization(context, context.config());
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // InitializeData

#endif // PRECICE_NO_MPI

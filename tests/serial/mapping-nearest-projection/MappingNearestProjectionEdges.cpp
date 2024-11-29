#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestProjection)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_DATA_TEST_CASE(MappingNearestProjectionEdges,
                     boost::unit_test::data::make({true, false}) * boost::unit_test::data::make({true, false}),
                     defineEdgesExplicitly, useBulkFunctions)
{
  PRECICE_TEST();
  /**
   * @brief Tests the Nearest Projection Mapping between two participants with explicit definition of edges
   *
   */
  testMappingNearestProjection(defineEdgesExplicitly, useBulkFunctions, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingNearestProjection

#endif // PRECICE_NO_MPI

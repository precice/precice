#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/Participant.hpp>
#include "helpers.hpp"

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuad
 */
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestProjection)
BOOST_DATA_TEST_CASE(QuadMappingDiagonalNearestProjectionEdgesWideKite,
                     boost::unit_test::data::make({true, false}) * boost::unit_test::data::make({true, false}),
                     defineEdgesExplicitly, useBulkFunctions)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  testQuadMappingNearestProjectionWideKite(defineEdgesExplicitly, useBulkFunctions, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingNearestProjection

#endif // PRECICE_NO_MPI

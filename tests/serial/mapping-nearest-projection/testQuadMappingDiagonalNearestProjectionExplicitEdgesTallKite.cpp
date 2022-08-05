#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include "helpers.hpp"

/**
 * @brief Tests the Nearest Projection Mapping on a single participant on a quad mesh of a tall kite with setMeshQuadWithEdges
 */
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestProjection)
BOOST_AUTO_TEST_CASE(testQuadMappingDiagonalNearestProjectionExplicitEdgesTallKite)
{
  PRECICE_TEST("SolverOne"_on(1_rank));
  bool defineEdgesExplicitly = true;
  testQuadMappingNearestProjectionTallKite(defineEdgesExplicitly, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingNearestProjection

#endif // PRECICE_NO_MPI

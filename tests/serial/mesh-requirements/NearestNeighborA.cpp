#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MeshRequirements)
BOOST_AUTO_TEST_CASE(NearestNeighborA)
{
  PRECICE_TEST(1_rank);
  precice::SolverInterface interface("A", context.config(), 0, 1);
  auto                     meshID = interface.getMeshID("MeshA");
  BOOST_TEST(!interface.isMeshConnectivityRequired(meshID));
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MeshRequirements

#endif // PRECICE_NO_MPI

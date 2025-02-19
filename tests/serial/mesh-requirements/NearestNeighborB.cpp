#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MeshRequirements)
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NearestNeighborB)
{
  PRECICE_TEST();
  precice::Participant interface("B", context.config(), 0, 1);
  auto                 meshName = "MeshB";
  BOOST_TEST(!interface.requiresMeshConnectivityFor(meshName));
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MeshRequirements

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(StationaryMappingWithSolverMesh)
BOOST_AUTO_TEST_CASE(StationaryMappingWithSolverMesh2D)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank));
  runTestStationaryMappingWithSolverMesh(context.config(), 2, context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // StationaryMappingWithSolvermesh

#endif // PRECICE_NO_MPI

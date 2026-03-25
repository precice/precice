#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Lifecycle)
// Test that triggering an error in initialize() (by not providing required initial data)
// sets the erroneous state and prevents further API usage.
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ErrorInInitialize)
{
  PRECICE_TEST();
  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshName = "MeshOne";
    double coords[] = {0.1, 1.2, 2.3};
    interface.setMeshVertex(meshName, coords);
    // intentionally skip writing initial data to trigger an error.
  } else {
    auto   meshName = "MeshTwo";
    double coords[] = {0.12, 1.21, 2.2};
    interface.setMeshVertex(meshName, coords);
  }

  // initialize() should throw because initial data was not provided
  BOOST_CHECK_THROW(interface.initialize(), ::precice::Error);

  // After the error, further API calls should be blocked
  BOOST_CHECK_THROW(interface.isCouplingOngoing(), ::precice::Error);

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Lifecycle
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

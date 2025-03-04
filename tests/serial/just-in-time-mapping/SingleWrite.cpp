#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <numeric>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a just-in-time mapping on one participant to a mesh defined
// by another participant.
// Here, we test that a just-in-time mapping must not conflict with a mapping,
// i.e., we define a mapping from a provided mesh to a received mesh
// and want to use a just-in-time write- mapping on the received mesh using the
// same data field: this would be conflicting, since there is no way to read
// either of the data fields on any other participant
BOOST_AUTO_TEST_CASE(SingleWrite)
{
  PRECICE_TEST();
  // Set up Participant
  BOOST_CHECK_THROW(precice::Participant couplingInterface(context.name, context.config(), 0, 1), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

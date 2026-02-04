#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(FirstParticipant)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(InfiniteAdvance)
{
  PRECICE_TEST();

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<double> coords{0.0, 0.0, 0.0};
  participant.setMeshVertex(context.name + "-Mesh", coords);

  if (context.isNamed("SolverOne")) {
    participant.initialize();

    // fails as max dt is infinite / max double
    BOOST_CHECK_EXCEPTION(participant.advance(participant.getMaxTimeStepSize()),
                          ::precice::Error,
                          ::precice::testing::errorContains("time-window-size method="));
  } else {
    // fails as the connection is cut
    BOOST_CHECK_THROW(participant.initialize(), ::precice::Error);
  }
}

BOOST_AUTO_TEST_SUITE_END() // FirstParticipant
BOOST_AUTO_TEST_SUITE_END() // Fundamental
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

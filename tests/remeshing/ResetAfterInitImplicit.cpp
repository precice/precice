#ifndef PRECICE_NO_MPI

#include "precice/Participant.hpp"
#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(ResetAfterInitImplicit)
{
  PRECICE_TEST();
  using namespace precice::testing;
  constexpr double     y = 0.0;
  precice::Participant participant{context.name, context.config(), context.rank, context.size};

  auto qt = [&] {
    if (context.isNamed("A")) {
      return QuickTest(participant, "MA"_mesh, "DB"_read, "DA"_write);
    }
    return QuickTest(participant, "MB"_mesh, "DA"_read, "DB"_write);
  }();

  qt.setVertices({0.0, y, 1.0, y})
      .initialize()
      // Reset mesh directly after initialization (issue #2093)
      .resetMesh()
      .setVertices({0.0, y, 1.0, y})
      .writeCheckpoint()
      .readCheckpoint()
      .advance()
      // TW 1 It 2
      .writeCheckpoint()
      .expect({0.00, 0.00})
      .readCheckpoint()
      .advance()
      // TW 2 It 1 time window complete
      .writeCheckpoint()
      .readCheckpoint()
      .advance()
      // TW 2 It 2
      .writeCheckpoint()
      .expect({0.00, 0.00})
      .readCheckpoint()
      .advance()
      // Done
      .expectCouplingCompleted();
  qt.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

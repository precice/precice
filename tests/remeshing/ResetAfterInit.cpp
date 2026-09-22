#ifndef PRECICE_NO_MPI

#include "precice/Participant.hpp"
#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(ResetAfterInit)
{
  PRECICE_TEST();
  using namespace precice::testing;
  constexpr double     y = 0.0;
  precice::Participant participant{context.name, context.config(), context.rank, context.size};

  // A - Writer, remeshes directly after initialization
  if (context.isNamed("A")) {
    QuickTest(participant, "MA"_mesh, "D"_write)
        .setVertices({0.0, y, 1.0, y})
        .initialize()
        // Reset mesh directly after initialization (issue #2093)
        .resetMesh()
        .setVertices({0.0, y, 1.0, y, 2.0, y})
        .write({0.01, 0.02, 0.03})
        .advance()
        .write({0.11, 0.12, 0.13})
        .advance()
        .finalize();
  }
  // B - Reader, static geometry
  if (context.isNamed("B")) {
    QuickTest(participant, "MB"_mesh, "D"_read)
        .setVertices({0.0, y, 1.0, y, 2.0, y})
        .initialize()
        .advance()
        .expect({0.01, 0.02, 0.03})
        .advance()
        .expect({0.11, 0.12, 0.13})
        .finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

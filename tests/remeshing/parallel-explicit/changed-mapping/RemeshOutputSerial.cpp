#ifndef PRECICE_NO_MPI

#include "precice/Participant.hpp"
#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelExplicit)
BOOST_AUTO_TEST_SUITE(ChangedMapping)
BOOST_AUTO_TEST_CASE(RemeshOutputSerial)
{
  using namespace precice::testing;
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));
  constexpr double y = 0.0;

  precice::Participant participant{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    QuickTest(participant, "MA"_mesh, "D"_write)
        .setVertices({0.0, y, 1.0, y})
        .initialize()
        .write({0.01, 0.02})
        .advance()
        .write({0.11, 0.12})
        .advance()
        .finalize();
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    QuickTest(participant, "MB"_mesh, "D"_read)
        .setVertices({0.0, y, 1.0, y})
        .initialize()
        .advance()
        .expect({0.01, 0.02})
        .resetMesh()
        .setVertices({0.0, y, 1.0, y, 2.0, y})
        .advance()
        .expect({0.11, 0.12, 0.12})
        .finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // ChangedMapping
BOOST_AUTO_TEST_SUITE_END() // ParallelExplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

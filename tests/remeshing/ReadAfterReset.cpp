#ifndef PRECICE_NO_MPI

#include "precice/Participant.hpp"
#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadAfterReset)
{
  PRECICE_TEST();
  using namespace precice::testing;
  constexpr double     y = 0.0;
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
    auto qt = QuickTest(participant, "MB"_mesh, "D"_read)
                  .setVertices({0.0, y, 1.0, y})
                  .initialize()
                  .advance()
                  .expect({0.01, 0.02})

                  .resetMesh()
                  .setVertices({0.0, y});

    BOOST_CHECK_THROW(qt.read(), ::precice::Error);

    qt.advance()
        .expect({0.11})
        .finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

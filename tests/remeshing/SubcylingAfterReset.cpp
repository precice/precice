#ifndef PRECICE_NO_MPI

#include "precice/Participant.hpp"
#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_CASE(SubcylingAfterReset)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));
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
      // TW 1 It 1 currently not allowed
      .writeCheckpoint();
  qt.readCheckpoint()
      .advance()
      // TW 1 It 2
      .writeCheckpoint()
      .expect({0.00, 0.00});
  qt.readCheckpoint()
      .advance()
      // TW 2 It 1 time window complete
      .writeCheckpoint()
      .resetMesh()
      .setVertices({0.0, y, 1.0, y})
      .readCheckpoint();
  BOOST_CHECK_THROW(qt.advance(0.5), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelExplicit)
BOOST_AUTO_TEST_SUITE(TwoMeshes)
BOOST_AUTO_TEST_CASE(LocalMapping)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  std::vector<double>   coords{1.0, 0.0, 2.0, 0.0}, coords2{0.5, 0.0, 1.5, 0.0, 2.5, 0.0};
  std::array<double, 2> data, expected;
  std::array<double, 3> data2, expected2;

  if (context.isNamed("A")) {
    std::vector<precice::VertexID> a1ids(2), a2ids(2);
    p.setMeshVertices("MA1", coords, a1ids);
    p.setMeshVertices("MA2", coords, a2ids);
    p.initialize();

    data.fill(-1);
    expected.fill(0);
    p.readData("MA1", "DB", a1ids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(0);
    p.readData("MA2", "DB", a2ids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    data.fill(1.1);
    p.writeData("MA1", "DA1", a1ids, data);
    data.fill(1.2);
    p.writeData("MA2", "DA2", a2ids, data);

    p.advance(p.getMaxTimeStepSize());

    // We remesh MA1 in this time window

    data.fill(-1);
    expected.fill(1);
    p.readData("MA1", "DB", a1ids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(1);
    p.readData("MA2", "DB", a2ids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    // Remesh MA1
    p.resetMesh("MA1");
    a1ids.resize(3);
    p.setMeshVertices("MA1", coords2, a1ids);

    data2.fill(2.1);
    p.writeData("MA1", "DA1", a1ids, data2);
    data.fill(2.2);
    p.writeData("MA2", "DA2", a2ids, data);

    p.advance(p.getMaxTimeStepSize());

    // Normal time window
    data2.fill(-1);
    expected2.fill(2);
    p.readData("MA1", "DB", a1ids, p.getMaxTimeStepSize(), data2);
    BOOST_TEST(data2 == expected2, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(2);
    p.readData("MA2", "DB", a2ids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    data2.fill(3.1);
    p.writeData("MA1", "DA1", a1ids, data2);
    data.fill(3.2);
    p.writeData("MA2", "DA2", a2ids, data);

    p.advance(p.getMaxTimeStepSize());

    BOOST_REQUIRE(!p.isCouplingOngoing());

    data2.fill(-1);
    expected2.fill(3);
    p.readData("MA1", "DB", a1ids, 0, data2);
    BOOST_TEST(data2 == expected2, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(3);
    p.readData("MA2", "DB", a2ids, 0, data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

  } else {
    std::vector<precice::VertexID> bids(2);
    p.setMeshVertices("MB", coords, bids);
    p.initialize();

    data.fill(-1);
    expected.fill(0);
    p.readData("MB", "DA1", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    p.readData("MB", "DA2", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    data.fill(1);
    p.writeData("MB", "DB", bids, data);

    p.advance(p.getMaxTimeStepSize());

    data.fill(-1);
    expected.fill(1.1);
    p.readData("MB", "DA1", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(1.2);
    p.readData("MB", "DA2", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    data.fill(2);
    p.writeData("MB", "DB", bids, data);

    p.advance(p.getMaxTimeStepSize());

    data.fill(-1);
    expected.fill(2.1);
    p.readData("MB", "DA1", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(2.2);
    p.readData("MB", "DA2", bids, p.getMaxTimeStepSize(), data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    data.fill(3);
    p.writeData("MB", "DB", bids, data);

    p.advance(p.getMaxTimeStepSize());

    BOOST_REQUIRE(!p.isCouplingOngoing());

    data.fill(-1);
    expected.fill(3.1);
    p.readData("MB", "DA1", bids, 0, data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    data.fill(-1);
    expected.fill(3.2);
    p.readData("MB", "DA2", bids, 0, data);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_SUITE_END() // TwoMeshes
BOOST_AUTO_TEST_SUITE_END() // ParallelExplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(InitialData)
PRECICE_TEST_SETUP("PA"_on(1_rank), "PB"_on(1_rank))
BOOST_AUTO_TEST_CASE(InterleavedCreationWithGradients)
{
  PRECICE_TEST();

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("PA")) {
    BOOST_TEST(p.requiresInitialData());
    BOOST_TEST(p.requiresGradientDataFor("MA", "DA"));
    std::array<double, 1> d1{1}, d2{2}, d3{3}, d4{4};
    std::array<double, 2> g1{1, 1}, g2{2, 2}, g3{3, 3}, g4{4, 4};
    std::array<double, 2> p1{1.0, 0.0}, p2{2.0, 0.0}, p3{3.0, 0.0}, p4{4.0, 0.0};
    auto                  v1 = p.setMeshVertex("MA", p1);
    p.writeData("MA", "DA", {&v1, 1}, d1);
    p.writeGradientData("MA", "DA", {&v1, 1}, g1);
    auto v2 = p.setMeshVertex("MA", p2);
    p.writeData("MA", "DA", {&v2, 1}, d2);
    p.writeGradientData("MA", "DA", {&v2, 1}, g2);
    auto v3 = p.setMeshVertex("MA", p3);
    p.writeData("MA", "DA", {&v3, 1}, d3);
    p.writeGradientData("MA", "DA", {&v3, 1}, g3);
    auto v4 = p.setMeshVertex("MA", p4);
    p.writeData("MA", "DA", {&v4, 1}, d4);
    p.writeGradientData("MA", "DA", {&v4, 1}, g4);
    p.initialize();
  } else {
    std::array<precice::VertexID, 4> vids;
    std::vector<double>              pos{1.0, 1.0, 2.0, 1.0, 3.0, 1.0, 4.0, 1.0};
    p.setMeshVertices("MB", pos, vids);
    p.initialize();
    std::array<double, 4> data;
    p.readData("MB", "DA", vids, 0, data);
    std::vector<double> expected{2.0, 4.0, 6.0, 8.0};
    BOOST_TEST(data == expected, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_SUITE_END() // InitialData
BOOST_AUTO_TEST_SUITE_END() // Fundamental
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingVolume)

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(QuadRead2D)
{
  PRECICE_TEST();

  precice::Participant p(context.name, context.config(), 0, 1);

  if (context.isNamed("SolverOne")) {
    std::vector<precice::VertexID> ids(4);
    std::vector<double>            coords{
        0, 0, // 1
        0, 4, // 10
        4, 4, // 100
        4, 0  // 1000
    };
    p.setMeshVertices("MeshOne", coords, ids);
    p.setMeshQuads("MeshOne", ids);
    std::vector<double> data{1, 10, 100, 1000};
    p.initialize();
    p.writeData("MeshOne", "DataOne", ids, data);
    p.advance(p.getMaxTimeStepSize());
  } else {
    std::vector<precice::VertexID> ids(7);
    std::vector<double>            coords{
        0, 0,  // hit point 0,0
        5, 5,  // fallback to point 4,4
        5, -1, // fallback to point 4,0
        2, 2,  // middle of diagonal 0,0 and 4,4
        5, 2,  // middle outside right of quad
        1, 3,  // left triangle closest to 0,4
        3, 1,  // right triangle closest to 4,0
    };
    p.setMeshVertices("MeshTwo", coords, ids);
    p.initialize();
    std::vector<double> data(ids.size());
    p.readData("MeshTwo", "DataOne", ids, p.getMaxTimeStepSize(), data);
    p.advance(p.getMaxTimeStepSize());
    std::vector<double> expected{
        1, 100, 1000, // points
        50.5, 550,    // edges
        30.25, 525.25 // triangles
    };
    BOOST_TEST(data == expected, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_SUITE_END() // MappingVolume
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#include <Eigen/Core>
#include <sstream>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(TetrahedronTests)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BasicTetra)
{
  PRECICE_TEST();
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(0.0, 0.0, 1.0);

  Vertex v1(coords1, 0);
  Vertex v2(coords2, 1);
  Vertex v3(coords3, 2);
  Vertex v4(coords4, 3);

  Tetrahedron tetra(v1, v2, v3, v4);

  Vertex &v1ref = tetra.vertex(0);
  BOOST_TEST(v1ref.getID() == v1.getID());

  Vertex &v2ref = tetra.vertex(1);
  BOOST_TEST(v2ref.getID() == v2.getID());

  Vertex &v3ref = tetra.vertex(2);
  BOOST_TEST(v3ref.getID() == v3.getID());

  Vertex &v4ref = tetra.vertex(3);
  BOOST_TEST(v4ref.getID() == v4.getID());

  Vector3d center = tetra.getCenter();
  BOOST_TEST(testing::equals(center, (coords1 + coords2 + coords3 + coords4) / 4));

  // sqrt(11)/4
  constexpr double expectedRadius = 0.82915619758;
  BOOST_TEST(tetra.getEnclosingRadius() == expectedRadius);

  constexpr double expectedVolume = 1.0 / 6.0;
  BOOST_TEST(tetra.getVolume() == expectedVolume);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WeirdTetra)
{
  PRECICE_TEST();
  // Same as above, but with a vertex whose projection
  // on the opposing tetra is out ouf the tetrahedron.
  // Also, we give it a negative z-coordinate
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(-5.0, 10.0, -1.0);

  Vertex v1(coords1, 0);
  Vertex v2(coords2, 1);
  Vertex v3(coords3, 2);
  Vertex v4(coords4, 3);

  Tetrahedron tetra(v1, v2, v3, v4);

  Vertex &v1ref = tetra.vertex(0);
  BOOST_TEST(v1ref.getID() == v1.getID());

  Vertex &v2ref = tetra.vertex(1);
  BOOST_TEST(v2ref.getID() == v2.getID());

  Vertex &v3ref = tetra.vertex(2);
  BOOST_TEST(v3ref.getID() == v3.getID());

  Vertex &v4ref = tetra.vertex(3);
  BOOST_TEST(v4ref.getID() == v4.getID());

  Vector3d center = tetra.getCenter();
  BOOST_TEST(testing::equals(center, (coords1 + coords2 + coords3 + coords4) / 4));

  // sqrt(11)/4
  constexpr double expectedRadius = 8.314144574157945;
  BOOST_TEST(tetra.getEnclosingRadius() == expectedRadius);

  constexpr double expectedVolume = 1.0 / 6.0;
  BOOST_TEST(tetra.getVolume() == expectedVolume);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TetrahedronEquality)
{
  PRECICE_TEST();
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(0.0, 0.0, 1.0);
  Vector3d coords5(0.0, 0.0, -1.0);

  Vertex v1(coords1, 0);
  Vertex v2(coords2, 1);
  Vertex v3(coords3, 2);
  Vertex v4(coords4, 3);
  Vertex v5(coords5, 4);

  Tetrahedron tetra1(v1, v2, v3, v4);
  Tetrahedron tetra2(v3, v1, v2, v4);
  Tetrahedron tetra3(v1, v2, v3, v5);

  BOOST_TEST(tetra1 == tetra2);
  BOOST_TEST(tetra1 != tetra3);
  BOOST_TEST(tetra2 != tetra3);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TetrahedronWKTPrint)
{
  PRECICE_TEST();
  Vertex v1(Eigen::Vector3d(0., 0., 0.), 0);
  Vertex v2(Eigen::Vector3d(1., 0., 0.), 1);
  Vertex v3(Eigen::Vector3d(0., 1., 0.), 2);
  Vertex v4(Eigen::Vector3d(0., 0., 1.), 3);

  Tetrahedron       t1(v1, v2, v3, v4);
  std::stringstream stream;
  stream << t1;
  std::string t1string("MULTILINESTRING ((0 0 0, 1 0 0), (0 0 0, 0 1 0), (0 0 0, 0 0 1), (1 0 0, 0 1 0), (1 0 0, 0 0 1), (0 1 0, 0 0 1))");
  BOOST_TEST(t1string == stream.str());
}

BOOST_AUTO_TEST_SUITE_END() // Tetrahedron
BOOST_AUTO_TEST_SUITE_END() // Mesh

#include <Eigen/Core>
#include <iterator>
#include <sstream>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/RangeAccessor.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(TriangleTests)

BOOST_AUTO_TEST_CASE(BasicTetra)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(0.0, 0.0, 1.0);

  Vertex v1(coords1, 0);
  Vertex v2(coords2, 1);
  Vertex v3(coords3, 2);
  Vertex v4(coords4, 3);

  Tetrahedron tetra(v1, v2, v3, v4, 0);

  Vertex &v1ref = tetra.vertex(0);
  BOOST_TEST(v1ref.getID() == v1.getID());

  Vertex &v2ref = tetra.vertex(1);
  BOOST_TEST(v2ref.getID() == v2.getID());

  Vertex &v3ref = tetra.vertex(2);
  BOOST_TEST(v3ref.getID() == v3.getID());

  Vertex &v4ref = tetra.vertex(3);
  BOOST_TEST(v4ref.getID() == v4.getID());

  TetrahedronID id = tetra.getID();
  BOOST_TEST(id == 0);

  Vector3d center = tetra.getCenter();
  BOOST_TEST(testing::equals(center, (coords1 + coords2 + coords3 + coords4) / 4));

  // sqrt(11)/4
  constexpr double expectedRadius = 0.82915619758;
  BOOST_TEST(tetra.getEnclosingRadius() == expectedRadius);

  constexpr double expectedVolume = 1.0 / 6.0;
  BOOST_TEST(tetra.getVolume() == expectedVolume);

}

BOOST_AUTO_TEST_CASE(WeirdTetra)
{

  // Same as above, but with a vertex whose projection
  // on the opposing tetra is out ouf the tetrahedron.
  // Also, we give it a negative z-coordinate
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(-5.0, 10.0, -1.0);

  Vertex v1(coords1, 0);
  Vertex v2(coords2, 1);
  Vertex v3(coords3, 2);
  Vertex v4(coords4, 3);

  Tetrahedron tetra(v1, v2, v3, v4, 0);

  Vertex &v1ref = tetra.vertex(0);
  BOOST_TEST(v1ref.getID() == v1.getID());

  Vertex &v2ref = tetra.vertex(1);
  BOOST_TEST(v2ref.getID() == v2.getID());

  Vertex &v3ref = tetra.vertex(2);
  BOOST_TEST(v3ref.getID() == v3.getID());

  Vertex &v4ref = tetra.vertex(3);
  BOOST_TEST(v4ref.getID() == v4.getID());

  TetrahedronID id = tetra.getID();
  BOOST_TEST(id == 0);

  Vector3d center = tetra.getCenter();
  BOOST_TEST(testing::equals(center, (coords1 + coords2 + coords3 + coords4) / 4));

  // sqrt(11)/4
  constexpr double expectedRadius = 8.314144574157945;
  BOOST_TEST(tetra.getEnclosingRadius() == expectedRadius);

  constexpr double expectedVolume = 1.0 / 6.0;
  BOOST_TEST(tetra.getVolume() == expectedVolume);

}

BOOST_AUTO_TEST_CASE(TetraRangeAccess)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d coords1(0.0, 0.0, 0.0);
  Vector3d coords2(1.0, 0.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(0.0, 0.0, 1.0);

  Vertex v0(coords1, 0);
  Vertex v1(coords2, 1);
  Vertex v2(coords3, 2);
  Vertex v3(coords4, 3);

  Tetrahedron tetra(v0, v1, v2, v3, 0);

  {
    // Test begin(), end()
    auto       ibegin = tetra.begin();
    const auto iend   = tetra.end();
    BOOST_TEST(std::distance(ibegin, iend) == 4);
    BOOST_TEST(*ibegin == v0.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v1.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v2.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v3.rawCoords());
    ++ibegin;
    BOOST_TEST((ibegin == iend));
  }
  {
    // Test begin(), end() for const
    const Tetrahedron &ctetra = tetra;
    auto            ibegin    = ctetra.begin();
    const auto      iend      = ctetra.end();
    BOOST_TEST(std::distance(ibegin, iend) == 4);
    BOOST_TEST(*ibegin == v0.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v1.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v2.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v3.rawCoords());
    ++ibegin;
    BOOST_TEST((ibegin == iend));
  }
  {
    // Test cbegin(), cend()
    auto       ibegin = tetra.cbegin();
    const auto iend   = tetra.cend();
    BOOST_TEST(std::distance(ibegin, iend) == 4);
    BOOST_TEST(*ibegin == v0.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v1.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v2.rawCoords());
    ++ibegin;
    BOOST_TEST(*ibegin == v3.rawCoords());
    ++ibegin;
    BOOST_TEST((ibegin == iend));
  }
}

BOOST_AUTO_TEST_CASE(TetrahedronEquality)
{
  PRECICE_TEST(1_rank);
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

  Tetrahedron tetra1(v1, v2, v3, v4, 0);
  Tetrahedron tetra2(v3, v1, v2, v4, 0);
  Tetrahedron tetra3(v1, v2, v3, v5, 0);

  BOOST_TEST(tetra1 == tetra2);
  BOOST_TEST(tetra1 != tetra3);
  BOOST_TEST(tetra2 != tetra3);
}
/*
BOOST_AUTO_TEST_CASE(TriangleWKTPrint)
{
  PRECICE_TEST(1_rank);
  Vertex            v1(Eigen::Vector3d(0., 0., 0.), 0);
  Vertex            v2(Eigen::Vector3d(0., 1., 0.), 0);
  Vertex            v3(Eigen::Vector3d(1., 0., 0.), 0);
  Edge              e1(v1, v2, 0);
  Edge              e2(v2, v3, 0);
  Edge              e3(v3, v1, 0);
  Triangle          t1(e1, e2, e3, 0);
  std::stringstream stream;
  stream << t1;
  std::string t1string("POLYGON ((0 0 0, 0 1 0, 1 0 0, 0 0 0))");
  BOOST_TEST(t1string == stream.str());
}*/

BOOST_AUTO_TEST_SUITE_END() // Triangle
BOOST_AUTO_TEST_SUITE_END() // Mesh

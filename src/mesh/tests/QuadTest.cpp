#include <sstream>

#include "mesh/Edge.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

#include <vector>

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(QuadTests)

BOOST_AUTO_TEST_CASE(Quads)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d coords0(0.0, 0.0, 0.0);
  Vector3d coords1(1.0, 0.0, 0.0);
  Vector3d coords2(1.0, 1.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v1, v2, 1);
    Edge e2(v2, v3, 2);
    Edge e3(v3, v0, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex &v0ref = quad.vertex(0);
    BOOST_TEST(v0ref.getID() == v0.getID());
    Vertex &v1ref = quad.vertex(1);
    BOOST_TEST(v1ref.getID() == v1.getID());
    Vertex &v2ref = quad.vertex(2);
    BOOST_TEST(v2ref.getID() == v2.getID());
    Vertex &v3ref = quad.vertex(3);
    BOOST_TEST(v3ref.getID() == v3.getID());

    Edge &e0ref = quad.edge(0);
    BOOST_TEST(e0ref.getID() == e0.getID());
    Edge &e1ref = quad.edge(1);
    BOOST_TEST(e1ref.getID() == e1.getID());
    Edge &e2ref = quad.edge(2);
    BOOST_TEST(e2ref.getID() == e2.getID());
    Edge &e3ref = quad.edge(3);
    BOOST_TEST(e3ref.getID() == e3.getID());

    int id = quad.getID();
    BOOST_TEST(id == 0);
  }
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v1, v2, 1);
    Edge e2(v3, v2, 2);
    Edge e3(v0, v3, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex &v0ref = quad.vertex(0);
    BOOST_TEST(v0ref.getID() == v0.getID());
    Vertex &v1ref = quad.vertex(1);
    BOOST_TEST(v1ref.getID() == v1.getID());
    Vertex &v2ref = quad.vertex(2);
    BOOST_TEST(v2ref.getID() == v2.getID());
    Vertex &v3ref = quad.vertex(3);
    BOOST_TEST(v3ref.getID() == v3.getID());

    Edge &e0ref = quad.edge(0);
    BOOST_TEST(e0ref.getID() == e0.getID());
    Edge &e1ref = quad.edge(1);
    BOOST_TEST(e1ref.getID() == e1.getID());
    Edge &e2ref = quad.edge(2);
    BOOST_TEST(e2ref.getID() == e2.getID());
    Edge &e3ref = quad.edge(3);
    BOOST_TEST(e3ref.getID() == e3.getID());

    int id = quad.getID();
    BOOST_TEST(id == 0);
  }
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v2, v1, 1);
    Edge e2(v2, v3, 2);
    Edge e3(v0, v3, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex &v0ref = quad.vertex(0);
    BOOST_TEST(v0ref.getID() == v0.getID());
    Vertex &v1ref = quad.vertex(1);
    BOOST_TEST(v1ref.getID() == v1.getID());
    Vertex &v2ref = quad.vertex(2);
    BOOST_TEST(v2ref.getID() == v2.getID());
    Vertex &v3ref = quad.vertex(3);
    BOOST_TEST(v3ref.getID() == v3.getID());

    Edge &e0ref = quad.edge(0);
    BOOST_TEST(e0ref.getID() == e0.getID());
    Edge &e1ref = quad.edge(1);
    BOOST_TEST(e1ref.getID() == e1.getID());
    Edge &e2ref = quad.edge(2);
    BOOST_TEST(e2ref.getID() == e2.getID());
    Edge &e3ref = quad.edge(3);
    BOOST_TEST(e3ref.getID() == e3.getID());

    int id = quad.getID();
    BOOST_TEST(id == 0);
  }
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v1, v2, 1);
    Edge e2(v2, v3, 2);
    Edge e3(v3, v0, 3);

    Quad quad(e0, e1, e2, e3, 0);
    {
      // Test begin(), end()
      auto       ibegin = quad.begin();
      const auto iend   = quad.end();
      BOOST_TEST(std::distance(ibegin, iend) == 4);
      BOOST_TEST(*ibegin == v0.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v1.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v2.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v3.getCoords());
      ++ibegin;
      BOOST_TEST((ibegin == iend));
    }
    {
      // Test begin(), end() for const
      const Quad &cquad  = quad;
      auto        ibegin = cquad.begin();
      const auto  iend   = cquad.end();
      BOOST_TEST(std::distance(ibegin, iend) == 4);
      BOOST_TEST(*ibegin == v0.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v1.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v2.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v3.getCoords());
      ++ibegin;
      BOOST_TEST((ibegin == iend));
    }
    {
      // Test cbegin(), cend()
      auto       ibegin = quad.cbegin();
      const auto iend   = quad.cend();
      BOOST_TEST(std::distance(ibegin, iend) == 4);
      BOOST_TEST(*ibegin == v0.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v1.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v2.getCoords());
      ++ibegin;
      BOOST_TEST(*ibegin == v3.getCoords());
      ++ibegin;
      BOOST_TEST((ibegin == iend));
    }
  }
}
BOOST_AUTO_TEST_CASE(QuadEquality)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d coords0(0.0, 0.0, 0.0);
  Vector3d coords1(1.0, 0.0, 0.0);
  Vector3d coords2(1.0, 1.0, 0.0);
  Vector3d coords3(0.0, 1.0, 0.0);
  Vector3d coords4(2.0, 0.0, 0.0);
  Vector3d coords5(2.0, 1.0, 0.0);

  Vertex v0(coords0, 0);
  Vertex v1(coords1, 1);
  Vertex v2(coords2, 2);
  Vertex v3(coords3, 3);
  Vertex v4(coords4, 0);
  Vertex v5(coords5, 1);

  Edge e0(v0, v1, 0);
  Edge e1(v2, v1, 1);
  Edge e2(v2, v3, 2);
  Edge e3(v0, v3, 3);
  Edge e4(v0, v4, 0);
  Edge e5(v4, v5, 1);
  Edge e6(v5, v3, 2);
  Edge e7(v3, v0, 3);

  //*  *
  //
  //*  *
  Quad quad1(e0, e1, e2, e3, 0);
  Quad quad1rev(e3, e2, e1, e0, 0);
  BOOST_TEST(quad1 == quad1rev);

  //*    *
  //
  //*    *
  Quad quad2(e4, e5, e6, e7, 0);
  Quad quad2n(e4, e5, e6, e7, 0);
  quad2n.setNormal(Vector3d(0., 0., 1.));
  BOOST_TEST(quad2 != quad1);
  BOOST_TEST(quad2 != quad2n);
}
BOOST_AUTO_TEST_CASE(QuadWKTPrint)
{
  PRECICE_TEST(1_rank);
  Vertex            v1(Eigen::Vector3d(0., 0., 0.), 0);
  Vertex            v2(Eigen::Vector3d(0., 1., 0.), 0);
  Vertex            v3(Eigen::Vector3d(1., 1., 0.), 0);
  Vertex            v4(Eigen::Vector3d(1., 0., 0.), 0);
  Edge              e1(v1, v2, 0);
  Edge              e2(v2, v3, 0);
  Edge              e3(v3, v4, 0);
  Edge              e4(v4, v1, 0);
  Quad              q1(e1, e2, e3, e4, 0);
  std::stringstream stream;
  stream << q1;
  std::string q1string("POLYGON ((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0))");
  BOOST_TEST(q1string == stream.str());
}

BOOST_AUTO_TEST_SUITE_END() // Quad
BOOST_AUTO_TEST_SUITE_END() // Mesh

#include "mesh/Edge.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

#include <vector>

using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(Quads)
{
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
      const Quad &cquad = quad;
      auto            ibegin    = cquad.begin();
      const auto      iend      = cquad.end();
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

BOOST_AUTO_TEST_SUITE_END() // Mesh

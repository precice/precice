#include "mesh/Utils.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(LocateInvalidId)
{
  using namespace precice;
  using namespace precice::mesh;
  mesh::Mesh mesh("2D Testmesh", 2, testing::nextMeshID());
  auto       v1 = mesh.createVertex(Eigen::Vector2d::Zero()).getID();
  auto       v2 = mesh.createVertex(Eigen::Vector2d::Ones()).getID();

  using VIDs = std::vector<VertexID>;

  {
    auto index = locateInvalidVertexID(mesh, VIDs{v1, v2, 99});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 2);
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{99, v1, v2});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 0);
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{99, v1, v2, 99});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 0);
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{v1, v2, 99, v1, v2});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 2);
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{99});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 0);
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{v1, v2});
    BOOST_TEST(!index.has_value());
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{});
    BOOST_TEST(!index.has_value());
  }
  {
    auto index = locateInvalidVertexID(mesh, VIDs{99});
    BOOST_TEST(index.has_value());
    BOOST_TEST(index.value() == 0);
  }
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE_END();

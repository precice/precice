#include <limits>
#include "mesh/Utils.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_SUITE(UtilsTests)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(LocateInvalidId)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FindDuplicateVertexCoords)
{
  PRECICE_TEST();
  using namespace precice;
  using namespace precice::mesh;
  mesh::Mesh mesh("3D Testmesh", 3, testing::nextMeshID());
  BOOST_TEST(findDuplicateVertexCoords(mesh).empty());

  mesh.createVertex(Eigen::Vector3d{0, 0, 0});
  mesh.createVertex(Eigen::Vector3d{1, 2, 3});
  mesh.createVertex(Eigen::Vector3d{1, 0, 0});
  BOOST_TEST(findDuplicateVertexCoords(mesh).empty());

  // (1,2,3) three times and (0,0,0) twice, each must be reported exactly once
  mesh.createVertex(Eigen::Vector3d{1, 2, 3});
  mesh.createVertex(Eigen::Vector3d{0, 0, 0});
  mesh.createVertex(Eigen::Vector3d{1, 2, 3});

  auto duplicates = findDuplicateVertexCoords(mesh);
  BOOST_TEST_REQUIRE(duplicates.size() == 2);
  BOOST_TEST(testing::equals(duplicates[0], Eigen::Vector3d{0, 0, 0}));
  BOOST_TEST(testing::equals(duplicates[1], Eigen::Vector3d{1, 2, 3}));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FindDuplicateVertexCoords2D)
{
  PRECICE_TEST();
  using namespace precice;
  using namespace precice::mesh;
  mesh::Mesh mesh("2D Testmesh", 2, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d{0.5, 1});
  mesh.createVertex(Eigen::Vector2d{1, 0.5});
  mesh.createVertex(Eigen::Vector2d{0.5, 1});

  auto duplicates = findDuplicateVertexCoords(mesh);
  BOOST_TEST_REQUIRE(duplicates.size() == 1);
  BOOST_TEST(duplicates[0].size() == 2);
  BOOST_TEST(testing::equals(duplicates[0], Eigen::Vector2d{0.5, 1}));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FindDuplicateVertexCoordsNonFinite)
{
  PRECICE_TEST();
  using namespace precice;
  using namespace precice::mesh;
  mesh::Mesh   mesh("3D Testmesh", 3, testing::nextMeshID());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 20; ++i) {
    mesh.createVertex(Eigen::Vector3d{static_cast<double>(i % 5), nan, 0});
  }
  mesh.createVertex(Eigen::Vector3d{inf, 0, 0});
  mesh.createVertex(Eigen::Vector3d{inf, 0, 0});
  mesh.createVertex(Eigen::Vector3d{1, 2, 3});
  mesh.createVertex(Eigen::Vector3d{1, 2, 3});

  // Non-finite vertices are ignored, finite duplicates are still found
  auto duplicates = findDuplicateVertexCoords(mesh);
  BOOST_TEST_REQUIRE(duplicates.size() == 1);
  BOOST_TEST(testing::equals(duplicates[0], Eigen::Vector3d{1, 2, 3}));
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE_END();

#include <Eigen/Core>
#include <algorithm>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(BoundingBoxTests)

BOOST_AUTO_TEST_CASE(ExpandByBoundingBox)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb1({0.0, 1.0,
                     0.0, 1.0,
                     0.0, 1.0});
    BoundingBox bb2({-1.0, 0.5,
                     2.0, 3.5,
                     0.0, 4.0});
    bb1.expandBy(bb2);
    std::vector<double> compareData = {-1.0, 1.0,
                                       0.0, 3.5,
                                       0.0, 4.0};
    BOOST_TEST(bb1.dataVector() == compareData);
  }
  { // 2D
    BoundingBox bb1({0.0, 2.0,
                     -1.0, 1.0});
    BoundingBox bb2({-1.0, 0.5,
                     2.0, 3.5});
    bb1.expandBy(bb2);
    std::vector<double> compareData = {-1.0, 2.0,
                                       -1.0, 3.5};
    BOOST_TEST(bb1.dataVector() == compareData);
  }
} // ExpandByBoundingBox

BOOST_AUTO_TEST_CASE(ExpandByVertex)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    0.0, 1.0,
                    0.0, 1.0});
    Vertex      v1(Eigen::Vector3d(-1.0, 3.0, 0.5), 0);
    bb.expandBy(v1);
    std::vector<double> compareData = {-1.0, 1.0,
                                       0.0, 3.0,
                                       0.0, 1.0};
    BOOST_TEST(bb.dataVector() == compareData);
  }
  { // 2D
    BoundingBox bb({-2.0, 1.0,
                    2.0, 4.0});
    Vertex      v1(Eigen::Vector2d(-4.0, 2.0), 0);
    bb.expandBy(v1);
    std::vector<double> compareData = {-4.0, 1.0,
                                       2.0, 4.0};
    BOOST_TEST(bb.dataVector() == compareData);
  }
} // ExpandByVertex

BOOST_AUTO_TEST_CASE(ExpandByRadius)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    0.5, 2.0,
                    -1.0, 3.0});
    double      supportRadius = 1.0;
    bb.expandBy(supportRadius);
    std::vector<double> compareData = {-1.0, 2.0,
                                       -0.5, 3.0,
                                       -2.0, 4.0};
    BOOST_TEST(bb.dataVector() == compareData);
  }
  { // 2D
    BoundingBox bb({-2.0, 1.0,
                    0.5, 2.0});
    double      supportRadius = 1.0;
    bb.expandBy(supportRadius);
    std::vector<double> compareData = {-3.0, 2.0,
                                       -0.5, 3.0};
    BOOST_TEST(bb.dataVector() == compareData);
  }
} // ExpandByRadius

BOOST_AUTO_TEST_CASE(Scaling)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb1({-1.0, 1.0,
                     0.5, 2.0,
                     1.0, 1.5});
    double      safetyFactor = 2.0;
    bb1.scaleBy(safetyFactor);
    std::vector<double> compareData = {-5.0, 5.0,
                                       -3.5, 6.0,
                                       -3.0, 5.5};
    BOOST_TEST(bb1.dataVector() == compareData);
  }
  { // 2D
    BoundingBox bb1({-1.0, 1.0,
                     0.5, 2.0});
    double      safetyFactor = 2.0;
    bb1.scaleBy(safetyFactor);
    std::vector<double> compareData = {-5.0, 5.0,
                                       -3.5, 6.0};
    BOOST_TEST(bb1.dataVector() == compareData);
  }
} // Scaling

BOOST_AUTO_TEST_CASE(CenterOfGravity)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0,
                    2.0, 4.0});

    Eigen::Vector3d compareCOG(0.5, 1.0, 3.0);
    BOOST_TEST(compareCOG == bb.center());
  }
  { // 2D
    BoundingBox bb({0.0, 1.0,
                    -2.0, 5.0});

    Eigen::Vector2d compareCOG(0.5, 1.5);
    BOOST_TEST(compareCOG == bb.center());
  }
} // CenterOfGravity

BOOST_AUTO_TEST_CASE(MinMaxCorner)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0,
                    2.0, 4.0});

    Eigen::Vector3d compareMin(0.0, -1.0, 2.0);
    Eigen::Vector3d compareMax(1.0, 3.0, 4.0);
    BOOST_TEST(compareMin == bb.minCorner());
    BOOST_TEST(compareMax == bb.maxCorner());
  }
  { // 2D
    BoundingBox bb({-1.0, 3.0,
                    2.0, 4.0});

    Eigen::Vector2d compareMin(-1.0, 2.0);
    Eigen::Vector2d compareMax(3.0, 4.0);
    BOOST_TEST(compareMin == bb.minCorner());
    BOOST_TEST(compareMax == bb.maxCorner());
  }
} // CenterOfGravity

BOOST_AUTO_TEST_CASE(Area)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0,
                    2.0, 4.0});
    {
      std::vector<bool> deadAxis    = {false, false, true};
      double            compareArea = 4.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
    {
      std::vector<bool> deadAxis    = {false, true, false};
      double            compareArea = 2.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
    {
      std::vector<bool> deadAxis    = {true, false, false};
      double            compareArea = 8.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
    {
      std::vector<bool> deadAxis    = {false, false, false};
      double            compareArea = 8.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
  }
  { // 2D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0});
    {
      std::vector<bool> deadAxis    = {false, true};
      double            compareArea = 1.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
    {
      std::vector<bool> deadAxis    = {true, false};
      double            compareArea = 4.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
    {
      std::vector<bool> deadAxis    = {false, false};
      double            compareArea = 4.0;
      BOOST_TEST(bb.getArea(deadAxis) == compareArea);
    }
  }
} // Area

BOOST_AUTO_TEST_CASE(Overlapping)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0,
                     2.0, 4.0});
    BoundingBox bb2({-1.0, 0.5,
                     2.0, 5.0,
                     1.0, 3.0});
    BoundingBox bb3({2.0, 5.0,
                     4.0, 5.0,
                     0.0, 1.0});

    BOOST_TEST(bb1.overlapping(bb2));
    BOOST_TEST(!bb1.overlapping(bb3));
  }
  { // 2D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0});
    BoundingBox bb2({-1.0, 0.5,
                     2.0, 5.0});
    BoundingBox bb3({2.0, 5.0,
                     4.0, 5.0});

    BOOST_TEST(bb1.overlapping(bb2));
    BOOST_TEST(!bb1.overlapping(bb3));
  }
} // Overalapping

BOOST_AUTO_TEST_CASE(Comparison)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0,
                     2.0, 4.0});
    BoundingBox bb2({0.0, 1.0,
                     -1.0, 3.0,
                     2.0, 4.0});
    BoundingBox bb3({2.0, 5.0,
                     4.0, 5.0,
                     0.0, 1.0});

    BOOST_TEST(bb1 == bb2);
    BOOST_TEST(!(bb1 == bb3));
  }
  { // 2D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0});
    BoundingBox bb2({0.0, 1.0,
                     -1.0, 3.0});
    BoundingBox bb3({2.0, 5.0,
                     4.0, 5.0});

    BOOST_TEST(bb1 == bb2);
    BOOST_TEST(!(bb1 == bb3));
  }
} // Comparison

BOOST_AUTO_TEST_CASE(Contains)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0,
                    2.0, 4.0});
    Vertex      v1(Eigen::Vector3d(0.2, 1.0, 3.0), 0);
    Vertex      v2(Eigen::Vector3d(1.2, -2.0, 5.0), 0);

    BOOST_TEST(bb.contains(v1));
    BOOST_TEST(!bb.contains(v2));
  }
  { // 2D
    BoundingBox bb({0.0, 1.0,
                    -1.0, 3.0});
    Vertex      v1(Eigen::Vector2d(0.2, 1.0), 0);
    Vertex      v2(Eigen::Vector2d(1.2, -2.0), 0);

    BOOST_TEST(bb.contains(v1));
    BOOST_TEST(!bb.contains(v2));
  }
  { // 3D Point
    BoundingBox bb({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    Vertex      v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
    Vertex      v2(Eigen::Vector3d(1.2, -2.0, 1.0), 0);

    BOOST_TEST(bb.contains(v1));
    BOOST_TEST(!bb.contains(v2));
  }
  { // 2D Point
    BoundingBox bb({0.0, 0.0, 0.0, 0.0});
    Vertex      v1(Eigen::Vector2d(0.0, 0.0), 0);
    Vertex      v2(Eigen::Vector2d(1.2, -2.0), 0);

    BOOST_TEST(bb.contains(v1));
    BOOST_TEST(!bb.contains(v2));
  }
} // Contains

BOOST_AUTO_TEST_CASE(EmptyCase)
{
  PRECICE_TEST(1_rank);
  { // 3D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0,
                     2.0, 4.0});
    BoundingBox bb2(3);

    BOOST_TEST(!bb1.empty());
    BOOST_TEST(bb2.empty());
  }
  { // 2D
    BoundingBox bb1({0.0, 1.0,
                     -1.0, 3.0});
    BoundingBox bb2(2);

    BOOST_TEST(!bb1.empty());
    BOOST_TEST(bb2.empty());
  }
} // EmptyCase

BOOST_AUTO_TEST_SUITE_END() // BoundingBox
BOOST_AUTO_TEST_SUITE_END() // Mesh

#include "mesh/BoundingBox.hpp"
#include "testing/Testing.hpp"
#include "utils/Helpers.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(BoundingBoxTests)

BOOST_AUTO_TEST_CASE(ExpandBy_BoundingBox)
{
  PRECICE_TEST(1_rank);
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

BOOST_AUTO_TEST_CASE(ExpandBy_Vertex)
{
  PRECICE_TEST(1_rank);
  BoundingBox bb({0.0, 1.0,
                   0.0, 1.0,
                   0.0, 1.0});
  Vertex v1(Eigen::Vector3d(-1.0, 3.0, 0.5), 0);
  bb.expandBy(v1);
  std::vector<double> compareData = {-1.0, 1.0,
                                      0.0, 3.0,
                                      0.0, 1.0};
  BOOST_TEST(bb.dataVector() == compareData);
}

BOOST_AUTO_TEST_CASE(ExpandBy_Radius)
{
  PRECICE_TEST(1_rank);
  BoundingBox bb({0.0, 1.0,
                   0.5, 2.0,
                  -1.0, 3.0});
  double supportRadius = 1.0;
  bb.expandBy(supportRadius);
  std::vector<double> compareData = {-1.0, 2.0,
                                     -0.5, 3.0,
                                     -2.0, 4.0};
  BOOST_TEST(bb.dataVector() == compareData);
}

BOOST_AUTO_TEST_CASE(Scaling)
{
  PRECICE_TEST(1_rank);
  BoundingBox bb1({-1.0, 1.0,
                   0.5, 2.0,
                   1.0, 1.5});
  double safetyFactor = 2.0;
  bb1.scaleBy(safetyFactor);
  std::vector<double> compareData = {-5.0, 5.0,
                                     -3.5, 6.0,
                                     -3.0, 5.5};
  BOOST_TEST(bb1.dataVector() == compareData);
}

BOOST_AUTO_TEST_CASE(CenterOfGravity)
{
  BoundingBox bb({0.0, 1.0,
                  -1.0, 3.0,
                   2.0, 4.0});

  Eigen::Vector3d compareCOG(0.5, 1.0, 3.0);
  BOOST_TEST(compareCOG == bb.center());
}

BOOST_AUTO_TEST_CASE(Area)
{
  BoundingBox bb({0.0, 1.0,
                 -1.0, 3.0,
                  2.0, 4.0});
  {
    std::vector<bool> deadAxis = {false, false, true};
    double compareArea = 4.0;
    BOOST_TEST(bb.getArea(deadAxis) == compareArea);
  }
  {
    std::vector<bool> deadAxis = {false, true, false};
    double compareArea = 2.0;
    BOOST_TEST(bb.getArea(deadAxis) == compareArea);
  }
  {
    std::vector<bool> deadAxis = {true, false, false};
    double compareArea = 8.0;
    BOOST_TEST(bb.getArea(deadAxis) == compareArea);
  }
}

BOOST_AUTO_TEST_CASE(Overlapping)
{
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

BOOST_AUTO_TEST_CASE(Comparison)
{
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

BOOST_AUTO_TEST_CASE(Contains)
{
  BoundingBox bb({0.0, 1.0,
                  -1.0, 3.0,
                   2.0, 4.0});
  Vertex v1(Eigen::Vector3d(0.2, 1.0, 3.0), 0);
  Vertex v2(Eigen::Vector3d(1.2, -2.0, 5.0), 0);

  BOOST_TEST(bb.contains(v1));
  BOOST_TEST(!bb.contains(v2));
}

BOOST_AUTO_TEST_CASE(EmptyCase)
{
  BoundingBox bb1({0.0, 1.0,
                  -1.0, 3.0,
                   2.0, 4.0});
  BoundingBox bb2(3);

  BOOST_TEST(!bb1.empty());
  BOOST_TEST(bb2.empty());
}

BOOST_AUTO_TEST_SUITE_END() // BoundingBox
BOOST_AUTO_TEST_SUITE_END() // Mesh
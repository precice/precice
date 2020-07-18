#include <Eigen/Core>
#include <ostream>
#include "math/barycenter.hpp"
#include "math/constants.hpp"
#include "math/differences.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::math::barycenter;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(Barycenter)

BOOST_AUTO_TEST_CASE(BarycenterEdge)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  using precice::testing::equals;
  Vector3d a(0.0, 0.0, 0.0);
  Vector3d b(1.0, 0.0, 0.0);
  Vector3d n(0.0, 0.0, 1.0);
  {
    Vector3d l(0.5, 0.0, 0.0);
    Vector2d coords(0.5, 0.5);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  {
    Vector3d l(0.0, 0.0, 0.0);
    Vector2d coords(1.0, 0.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  {
    Vector3d l(1.0, 0.0, 0.0);
    Vector2d coords(0, 1.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  {
    Vector3d l(0.75, 1.0, 0.0);
    Vector3d projected(0.75, 0.0, 0.0);
    Vector2d coords(0.25, 0.75);
    auto     ret = calcBarycentricCoordsForEdge(a, b, n, l);
    BOOST_TEST(equals(ret.projected, projected));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords), "Coords are " << ret.barycentricCoords << " but should be " << coords);
  }
}

BOOST_AUTO_TEST_CASE(BarycenterTriangle)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  using precice::testing::equals;
  Vector3d a(0.0, 0.0, 0.0);
  Vector3d b(0.0, 1.0, 0.0);
  Vector3d c(1.0, 0.0, 0.0);
  Vector3d n(0.0, 0.0, 1.0);
  // is a?
  {
    Vector3d coords(1.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, a);
    BOOST_TEST(equals(ret.projected, a));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is b?
  {
    Vector3d coords(0.0, 1.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, b);
    BOOST_TEST(equals(ret.projected, b));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is c?
  {
    Vector3d coords(0.0, 0.0, 1.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, c);
    BOOST_TEST(equals(ret.projected, c));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is in the middle
  {
    Vector3d l = (a + b + c) / 3;
    Vector3d coords(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is in the middle of ab
  {
    Vector3d l = (a + b) / 2;
    Vector3d coords(0.5, 0.5, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is in the middle of bc
  {
    Vector3d l = (b + c) / 2;
    Vector3d coords(0.0, 0.5, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is in the middle of ca
  {
    Vector3d l = (a + c) / 2;
    Vector3d coords(0.5, 0.0, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, l);
    BOOST_TEST(equals(ret.projected, l));
    BOOST_TEST(ret.barycentricCoords.sum() == 1.0);
    BOOST_TEST(equals(ret.barycentricCoords, coords));
  }
  // is outside
  {
    Vector3d l(2.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, n, l);
    BOOST_TEST((ret.barycentricCoords.array() < -precice::math::NUMERICAL_ZERO_DIFFERENCE).any(), "Min 1 coord should be negative " << ret.barycentricCoords);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Barycenter

BOOST_AUTO_TEST_SUITE_END() // Math

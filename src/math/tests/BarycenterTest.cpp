#include <Eigen/Core>
#include <boost/mpl/vector.hpp>
#include <ostream>
#include "math/barycenter.hpp"
#include "math/constants.hpp"
#include "math/differences.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::math::barycenter;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(BarycenterEdge)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BarycenterEdge2D)
{
  PRECICE_TEST();
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  using precice::testing::equals;
  Vector2d a(0.0, 0.0);
  Vector2d b(1.0, 0.0);
  Vector2d n(0.0, 1.0);
  {
    Vector2d l(0.5, 0.0);
    Vector2d coords(0.5, 0.5);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector2d l(0.0, 0.0);
    Vector2d coords(1.0, 0.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector2d l(1.0, 0.0);
    Vector2d coords(0, 1.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector2d l(0.75, 1.0);
    Vector2d coords(0.25, 0.75);
    auto     ret = calcBarycentricCoordsForEdge(a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords), "Coords are " << ret << " but should be " << coords);
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BarycenterEdge3D)
{
  PRECICE_TEST();
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
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector3d l(0.0, 0.0, 0.0);
    Vector2d coords(1.0, 0.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector3d l(1.0, 0.0, 0.0);
    Vector2d coords(0, 1.0);
    auto     ret = calcBarycentricCoordsForEdge(
        a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  {
    Vector3d l(0.75, 1.0, 0.0);
    Vector2d coords(0.25, 0.75);
    auto     ret = calcBarycentricCoordsForEdge(a, b, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords), fmt::format("Coords are {} but should be {}", ret, coords));
  }
}

BOOST_AUTO_TEST_SUITE_END() // BarycenterEdges

BOOST_AUTO_TEST_SUITE(BarycenterTriangle)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BarycenterTriangle3D)
{
  PRECICE_TEST();
  using Eigen::Vector3d;
  using precice::testing::equals;
  Vector3d a(0.0, 0.0, 0.0);
  Vector3d b(0.0, 1.0, 0.0);
  Vector3d c(1.0, 0.0, 0.0);
  Vector3d n(0.0, 0.0, 1.0);
  // is a?
  {
    Vector3d coords(1.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, a);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is b?
  {
    Vector3d coords(0.0, 1.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, b);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is c?
  {
    Vector3d coords(0.0, 0.0, 1.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, c);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle
  {
    Vector3d l = (a + b + c) / 3;
    Vector3d coords(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of ab
  {
    Vector3d l = (a + b) / 2;
    Vector3d coords(0.5, 0.5, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of bc
  {
    Vector3d l = (b + c) / 2;
    Vector3d coords(0.0, 0.5, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of ca
  {
    Vector3d l = (a + c) / 2;
    Vector3d coords(0.5, 0.0, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is outside
  {
    Vector3d l(2.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST((ret.array() < -precice::math::NUMERICAL_ZERO_DIFFERENCE).any(), fmt::format("Min 1 coord should be negative {}", ret));
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BarycenterTriangle2D)
{
  PRECICE_TEST();
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  using precice::testing::equals;
  Vector2d a(0.0, 0.0);
  Vector2d b(0.0, 1.0);
  Vector2d c(1.0, 0.0);
  Vector2d n(0.0, 0.0);
  // is a?
  {
    Vector3d coords(1.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, a);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is b?
  {
    Vector3d coords(0.0, 1.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, b);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is c?
  {
    Vector3d coords(0.0, 0.0, 1.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, c);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle
  {
    Vector2d l = (a + b + c) / 3;
    Vector3d coords(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of ab
  {
    Vector2d l = (a + b) / 2;
    Vector3d coords(0.5, 0.5, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of bc
  {
    Vector2d l = (b + c) / 2;
    Vector3d coords(0.0, 0.5, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is in the middle of ca
  {
    Vector2d l = (a + c) / 2;
    Vector3d coords(0.5, 0.0, 0.5);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // is outside
  {
    Vector2d l(2.0, 0.0);
    auto     ret = calcBarycentricCoordsForTriangle(a, b, c, l);
    BOOST_TEST((ret.array() < -precice::math::NUMERICAL_ZERO_DIFFERENCE).any(), fmt::format("Min 1 coord should be negative {}", ret));
  }
}

BOOST_AUTO_TEST_SUITE_END() // BarycenterTriangles

BOOST_AUTO_TEST_SUITE(BarycenterTetrahedra)

struct TetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, 1.0};

  Eigen::Vector3d center{0.25, 0.25, 0.25};
};

struct FlippedTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, -1.0};

  Eigen::Vector3d center{0.25, 0.25, -0.25};
};

struct AlmostDegenerateTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, 0.001};

  Eigen::Vector3d center{0.25, 0.25, 0.00025};
};

struct FunnyTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{11.0, 11.0, 1.0};

  Eigen::Vector3d center{3.0, 3.0, 0.25};
};

typedef boost::mpl::vector<TetrahedronFixture, FlippedTetrahedronFixture, AlmostDegenerateTetrahedronFixture, FunnyTetrahedronFixture> TetrahedraFixtures;

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExactOnA, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(1.0, 0.0, 0.0, 0.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::a);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExactOnB, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(0.0, 1.0, 0.0, 0.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::b);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExactOnC, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(0.0, 0.0, 1.0, 0.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::c);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExactOnD, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(0.0, 0.0, 0.0, 1.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::d);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExactOnCenter, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d center_coords{0.25, 0.25, 0.25, 0.25};
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::center);
  BOOST_TEST(precice::testing::equals(ret, center_coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronInsidePoint, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(0.2, 0.3, 0.4, 0.1);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, 0.2 * T::a + 0.3 * T::b + 0.4 * T::c + 0.1 * T::d);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronEdgeCenter, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();
  using Eigen::Vector4d;

  Eigen::Vector4d coords(0.5, 0.5, 0.0, 0.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, 0.5 * T::a + 0.5 * T::b);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronTriangleCenter, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();
  using Eigen::Vector4d;

  Eigen::Vector4d coords(1. / 3, 1. / 3, 0.0, 1. / 3);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, (T::a + T::b + T::d) / 3);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExtrapolationOnEdge, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(-2.0, 3.0, 0, 0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, T::a + 3 * (T::b - T::a));
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExtrapolationOnMirrored, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(2.0, 0.0, 0.0, -1.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, 2 * T::a - T::d);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExtrapolationOnMirroredFromTriangle, T, TetrahedraFixtures, T)
{
  PRECICE_TEST();

  Eigen::Vector4d coords(1.0, 1.0, 1.0, -2.0);
  auto            ret = calcBarycentricCoordsForTetrahedron(T::a, T::b, T::c, T::d, (T::a + T::b + T::c) - 2 * T::d);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(precice::testing::equals(ret, coords));
}

BOOST_AUTO_TEST_SUITE_END() // BarycenterTetrahedra

BOOST_AUTO_TEST_SUITE_END() // Math

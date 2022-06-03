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

BOOST_AUTO_TEST_CASE(BarycenterEdge2D)
{
  PRECICE_TEST(1_rank);
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

BOOST_AUTO_TEST_CASE(BarycenterEdge3D)
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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(BarycenterTriangle)

BOOST_AUTO_TEST_CASE(BarycenterTriangle3D)
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

BOOST_AUTO_TEST_CASE(BarycenterTriangle2D)
{
  PRECICE_TEST(1_rank);
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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(BarycenterTetrahedra)

struct TetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, 1.0};

  Eigen::Vector3d center{0.25, 0.25, 0.25};
  Eigen::Vector4d center_coords{0.25, 0.25, 0.25, 0.25};
};

struct FlippedTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, -1.0};

  Eigen::Vector3d center{0.25, 0.25, -0.25};
  Eigen::Vector4d center_coords{0.25, 0.25, 0.25, 0.25};
};

struct AlmostDegenerateTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{0.0, 0.0, 0.001};

  Eigen::Vector3d center{0.25, 0.25, 0.00025};
  Eigen::Vector4d center_coords{0.25, 0.25, 0.25, 0.25};
};

struct FunnyTetrahedronFixture {
  Eigen::Vector3d a{0.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{1.0, 0.0, 0.0};
  Eigen::Vector3d d{11.0, 11.0, 1.0};

  Eigen::Vector3d center{3.0, 3.0, 0.25};
};

typedef boost::mpl::vector<TetrahedronFixture, FlippedTetrahedronFixture, AlmostDegenerateTetrahedronFixture, FunnyTetrahedronFixture> TetrahedraFixtures;

BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronVertices, T, TetrahedraFixtures, T)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  using Eigen::Vector4d;
  using precice::testing::equals;

  // Alias for fixture data
  const auto &a = T::a;
  const auto &b = T::b;
  const auto &c = T::c;
  const auto &d = T::d;

  // Is A?
  {
    Vector4d coords(1.0, 0.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, a);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Is B?
  {
    Vector4d coords(0.0, 1.0, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, b);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Is C?
  {
    Vector4d coords(0.0, 0.0, 1.0, 0.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, c);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Is D?
  {
    Vector4d coords(0.0, 0.0, 0.0, 1.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, d);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronCenter, T, TetrahedraFixtures, T)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  using Eigen::Vector4d;
  using precice::testing::equals;

  // Alias for fixture data
  const auto &a      = T::a;
  const auto &b      = T::b;
  const auto &c      = T::c;
  const auto &d      = T::d;
  const auto &center = T::center;

  Vector4d center_coords{0.25, 0.25, 0.25, 0.25};
  auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, center);
  BOOST_TEST(ret.sum() == 1.0);
  BOOST_TEST(equals(ret, center_coords));
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronInterpolation, T, TetrahedraFixtures, T)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  using Eigen::Vector4d;
  using precice::testing::equals;

  // Alias for fixture data
  const auto &a = T::a;
  const auto &b = T::b;
  const auto &c = T::c;
  const auto &d = T::d;

  // Arbitrary combination
  {
    Vector4d coords(0.2, 0.3, 0.4, 0.1);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, 0.2 * a + 0.3 * b + 0.4 * c + 0.1 * d);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }

  // Middle of edge AB
  {
    Vector4d coords(0.5, 0.5, 0.0, 0.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, 0.5 * a + 0.5 * b);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Middle of triangle ABD
  {
    Vector4d coords(1. / 3, 1. / 3, 0.0, 1. / 3);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, (a + b + d) / 3);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(BarycenterTetrahedronExtrapolation, T, TetrahedraFixtures, T)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  using Eigen::Vector4d;
  using precice::testing::equals;

  // Alias for fixture data
  const auto &a = T::a;
  const auto &b = T::b;
  const auto &c = T::c;
  const auto &d = T::d;

  // A + 3*(B-A)
  {
    Vector4d coords(-2.0, 3.0, 0, 0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, a + 3 * (b - a));
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Mirror of D with respect to A
  {
    Vector4d coords(2.0, 0.0, 0.0, -1.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, 2 * a - d);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
  // Mirrored D with ABC as symmetry plane
  {
    Vector4d coords(1.0, 1.0, 1.0, -2.0);
    auto     ret = calcBarycentricCoordsForTetrahedron(a, b, c, d, (a + b + c) - 2 * d);
    BOOST_TEST(ret.sum() == 1.0);
    BOOST_TEST(equals(ret, coords));
  }
}

BOOST_AUTO_TEST_SUITE_END() // Barycenter

BOOST_AUTO_TEST_SUITE_END() // Math

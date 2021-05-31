#include <Eigen/Core>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/Dimensions.hpp"
#include "utils/String.hpp"

using namespace precice;
using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(DimensionTests)
BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(LinearizeDelinearize)
{
  PRECICE_TEST(1_rank);
  { // 2D
    using Eigen::Vector2d;
    BOOST_TEST(linearize(Vector2d(0.0, 0.0)) == 0);
    BOOST_TEST(linearize(Vector2d(1.0, 0.0)) == 1);
    BOOST_TEST(linearize(Vector2d(0.0, 1.0)) == 2);
    BOOST_TEST(linearize(Vector2d(1.0, 1.0)) == 3);

    Vector2d delin0 = delinearize(0, 2);
    BOOST_TEST(delin0(0) == -1.0);
    BOOST_TEST(delin0(1) == -1.0);

    Vector2d delin1 = delinearize(1, 2);
    BOOST_TEST(delin1(0) == 1.0);
    BOOST_TEST(delin1(1) == -1.0);

    Vector2d delin2 = delinearize(2, 2);
    BOOST_TEST(delin2(0) == -1.0);
    BOOST_TEST(delin2(1) == 1.0);

    Vector2d delin3 = delinearize(3, 2);
    BOOST_TEST(delin3(0) == 1.0);
    BOOST_TEST(delin3(1) == 1.0);
  }

  { // 3D
    using Eigen::Vector3d;
    BOOST_TEST(linearize(Vector3d(0.0, 0.0, 0.0)) == 0);
    BOOST_TEST(linearize(Vector3d(1.0, 0.0, 0.0)) == 1);
    BOOST_TEST(linearize(Vector3d(0.0, 1.0, 0.0)) == 2);
    BOOST_TEST(linearize(Vector3d(1.0, 1.0, 0.0)) == 3);
    BOOST_TEST(linearize(Vector3d(0.0, 0.0, 1.0)) == 4);
    BOOST_TEST(linearize(Vector3d(1.0, 0.0, 1.0)) == 5);
    BOOST_TEST(linearize(Vector3d(0.0, 1.0, 1.0)) == 6);
    BOOST_TEST(linearize(Vector3d(1.0, 1.0, 1.0)) == 7);

    Vector3d delin0 = delinearize(0, 3);
    BOOST_TEST(delin0(0) == -1.0);
    BOOST_TEST(delin0(1) == -1.0);
    BOOST_TEST(delin0(2) == -1.0);

    Vector3d delin1 = delinearize(1, 3);
    BOOST_TEST(delin1(0) == 1.0);
    BOOST_TEST(delin1(1) == -1.0);
    BOOST_TEST(delin1(2) == -1.0);

    Vector3d delin2 = delinearize(2, 3);
    BOOST_TEST(delin2(0) == -1.0);
    BOOST_TEST(delin2(1) == 1.0);
    BOOST_TEST(delin2(2) == -1.0);

    Vector3d delin3 = delinearize(3, 3);
    BOOST_TEST(delin3(0) == 1.0);
    BOOST_TEST(delin3(1) == 1.0);
    BOOST_TEST(delin3(2) == -1.0);

    Vector3d delin4 = delinearize(4, 3);
    BOOST_TEST(delin4(0) == -1.0);
    BOOST_TEST(delin4(1) == -1.0);
    BOOST_TEST(delin4(2) == 1.0);

    Vector3d delin5 = delinearize(5, 3);
    BOOST_TEST(delin5(0) == 1.0);
    BOOST_TEST(delin5(1) == -1.0);
    BOOST_TEST(delin5(2) == 1.0);

    Vector3d delin6 = delinearize(6, 3);
    BOOST_TEST(delin6(0) == -1.0);
    BOOST_TEST(delin6(1) == 1.0);
    BOOST_TEST(delin6(2) == 1.0);

    Vector3d delin7 = delinearize(7, 3);
    BOOST_TEST(delin7(0) == 1.0);
    BOOST_TEST(delin7(1) == 1.0);
    BOOST_TEST(delin7(2) == 1.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#include <Eigen/Core>
#include <algorithm>
#include <stddef.h>
#include <vector>
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/ResidualPreconditioner.hpp"
#include "acceleration/impl/ResidualSumPreconditioner.hpp"
#include "acceleration/impl/ValuePreconditioner.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using namespace precice;
using namespace precice::acceleration;
using namespace precice::acceleration::impl;

struct ResPreconditionerFixture {
  Eigen::VectorXd _data;
  Eigen::VectorXd _res;
  Eigen::VectorXd _compareDataRes;
  Eigen::VectorXd _compareDataResSum;
  Eigen::VectorXd _compareDataResSum2;
  Eigen::VectorXd _compareDataValue;
  Eigen::VectorXd _compareDataConstant;

  ResPreconditionerFixture()
  {
    _data.resize(8);
    _data << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;

    _res.resize(8);
    _res << 0.1, 0.1, 0.001, 0.001, 0.001, 0.001, 10.0, 20.0;

    _compareDataRes.resize(8);
    _compareDataRes << 7.07106781186547372897e+00,
        1.41421356237309474579e+01,
        1.50000000000000000000e+03,
        2.00000000000000000000e+03,
        2.50000000000000000000e+03,
        3.00000000000000000000e+03,
        3.13049516849970566046e-01,
        3.57770876399966353265e-01;

    _compareDataResSum.resize(8);
    _compareDataResSum << 7.90585229434499154877e+01,
        1.58117045886899830975e+02,
        1.67708453051717078779e+04,
        2.23611270735622783832e+04,
        2.79514088419528488885e+04,
        3.35416906103434157558e+04,
        3.50007001329973377324e+00,
        4.00008001519969536020e+00;

    _compareDataResSum2.resize(8);
    _compareDataResSum2 << 1.58113093108981217938e+02,
        3.16226186217962435876e+02,
        4.74339279326943596971e+02,
        4.00008000319945455914e+00,
        5.00010000399932064141e+00,
        6.00012000479918228280e+00,
        7.00014000559904481236e+00,
        8.00016000639890734192e+00;

    _compareDataValue.resize(8);
    _compareDataValue << 4.47213595499957927704e-01,
        8.94427190999915855407e-01,
        3.23498319610315276940e-01,
        4.31331092813753647075e-01,
        5.39163866017192239255e-01,
        6.46996639220630553879e-01,
        6.58504607868518165859e-01,
        7.52576694706877713514e-01;

    _compareDataConstant.resize(8);
    _compareDataConstant << 1.00000000000000002082e-03,
        2.00000000000000004163e-03,
        1.49999999999999977796e+00,
        1.99999999999999955591e+00,
        2.50000000000000044409e+00,
        2.99999999999999955591e+00,
        6.99999999999999883585e+05,
        7.99999999999999650754e+05;
  }
  ~ResPreconditionerFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(ResPreconditionerTests, ResPreconditionerFixture)

BOOST_AUTO_TEST_CASE(testResPreconditioner)
{
  PRECICE_TEST(1_rank);
  std::vector<size_t> svs;
  svs.push_back(2);
  svs.push_back(4);
  svs.push_back(2);

  ResidualPreconditioner precond(-1);

  precond.initialize(svs);
  Eigen::VectorXd backup = _data;

  //should change
  precond.update(false, _data, _res);
  BOOST_TEST(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataRes));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //should not change weights
  precond.update(true, _data, _res * 10);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataRes));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));
}

BOOST_AUTO_TEST_CASE(testResSumPreconditioner)
{
  PRECICE_TEST(1_rank);
  std::vector<size_t> svs;
  svs.push_back(2);
  svs.push_back(4);
  svs.push_back(2);

  ResidualSumPreconditioner precond(-1);

  precond.initialize(svs);
  Eigen::VectorXd backup = _data;

  //should change, update twice to really test the summation
  precond.update(false, _data, _res);
  precond.update(false, _data, _res * 2);
  BOOST_TEST(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataResSum));

  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //should not change weights
  precond.update(true, _data, _res * 10);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataResSum));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));
}

BOOST_AUTO_TEST_CASE(testValuePreconditioner)
{
  PRECICE_TEST(1_rank);
  std::vector<size_t> svs;
  svs.push_back(2);
  svs.push_back(4);
  svs.push_back(2);

  ValuePreconditioner precond(-1);

  precond.initialize(svs);
  Eigen::VectorXd backup = _data;

  //should change, since first timestep
  precond.update(false, _data, _res);
  BOOST_TEST(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataValue));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //now no change
  precond.update(false, _data, _res);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataValue));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //should change weights
  precond.update(true, _data * 2, _res);
  BOOST_TEST(precond.requireNewQR());
  precond.newQRfulfilled();
}

BOOST_AUTO_TEST_CASE(testConstPreconditioner)
{
  PRECICE_TEST(1_rank);
  std::vector<size_t> svs;
  svs.push_back(2);
  svs.push_back(4);
  svs.push_back(2);

  std::vector<double> factors;
  factors.push_back(1e3);
  factors.push_back(2.0);
  factors.push_back(1e-5);

  ConstantPreconditioner precond(factors);

  precond.initialize(svs); //new weights already computed here
  Eigen::VectorXd backup = _data;

  // should have no effect
  precond.update(false, _data, _res);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataConstant));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //should not change weights
  precond.update(true, _data, _res);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataConstant));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));
}

BOOST_AUTO_TEST_CASE(testMultilpleMeshes)
{
  PRECICE_TEST(1_rank);
  std::vector<size_t> svs;
  svs.push_back(3);
  svs.push_back(5);

  ResidualSumPreconditioner precond(-1);

  precond.initialize(svs);
  Eigen::VectorXd backup = _data;

  //should change
  precond.update(false, _data, _res);
  BOOST_TEST(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataResSum2));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));

  //should not change weights
  precond.update(true, _data, _res * 10);
  BOOST_TEST(not precond.requireNewQR());
  precond.apply(_data);
  BOOST_TEST(testing::equals(_data, _compareDataResSum2));
  precond.revert(_data);
  BOOST_TEST(testing::equals(_data, backup));
}

#ifndef PRECICE_NO_MPI
BOOST_AUTO_TEST_CASE(testParallelMatrixScaling)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  //setup data
  int localN = -1;
  if (context.isMaster()) {
    localN = 2;
  } else if (context.isRank(1)) {
    localN = 1;
  } else if (context.isRank(2)) {
    localN = 0;
  } else if (context.isRank(3)) {
    localN = 1;
  }

  int globalN = 4;

  Eigen::MatrixXd V(localN, 2);
  Eigen::MatrixXd M(globalN, localN);
  Eigen::VectorXd x(localN);
  Eigen::MatrixXd V_back(localN, 2);
  Eigen::MatrixXd M_back(globalN, localN);
  Eigen::VectorXd x_back(localN);

  if (context.isMaster()) {
    V(0, 0) = 1.0;
    V(0, 1) = 2.0;
    V(1, 0) = 3.0;
    V(1, 1) = 4.0;
    M(0, 0) = 1.0;
    M(0, 1) = 2.0;
    M(1, 0) = 3.0;
    M(1, 1) = 4.0;
    M(2, 0) = 1.0;
    M(2, 1) = 2.0;
    M(3, 0) = 3.0;
    M(3, 1) = 4.0;
    x(0)    = 5.0;
    x(1)    = 5.0;
  } else if (context.isRank(1)) {
    V(0, 0) = 5.0;
    V(0, 1) = 6.0;
    M(0, 0) = 1.0;
    M(1, 0) = 2.0;
    M(2, 0) = 3.0;
    M(3, 0) = 4.0;
    x(0)    = 5.0;
  } else if (context.isRank(2)) {
  } else if (context.isRank(3)) {
    V(0, 0) = 7.0;
    V(0, 1) = 8.0;
    M(0, 0) = 1.0;
    M(1, 0) = 2.0;
    M(2, 0) = 3.0;
    M(3, 0) = 4.0;
    x(0)    = 5.0;
  }

  V_back = V;
  M_back = M;
  x_back = x;

  std::vector<size_t> svs;
  svs.push_back(localN);

  ValuePreconditioner precond(-1);
  precond.initialize(svs);
  precond.update(true, x, x);
  BOOST_TEST(precond.requireNewQR());

  precond.apply(V);

  BOOST_TEST(testing::equals(V, V_back * 0.1));

  precond.revert(V);

  BOOST_TEST(testing::equals(V, V_back));
}
#endif

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

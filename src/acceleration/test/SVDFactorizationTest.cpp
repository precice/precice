#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>

#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "acceleration/impl/SVDFactorization.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using namespace precice;
using namespace precice::acceleration;
using namespace precice::acceleration::impl;

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testSVDFactorization)
{
  PRECICE_TEST();
  int                      m                   = 2; // size of a*
  int                      n                   = 3; // size of b*
  double                   eps                 = 1.2e-3;
  const bool               cyclicCommunication = false;
  std::vector<double>      factors(1.0, 1.0);
  Eigen::VectorXd          a1(2);
  Eigen::VectorXd          a2(2);
  Eigen::VectorXd          a3(2);
  Eigen::VectorXd          b1(3);
  Eigen::VectorXd          b2(3);
  Eigen::VectorXd          b3(3);
  SVDFactorization::Vector singleValues = SVDFactorization::Vector::Zero(2);

  a1 << 1., -1.;
  a2 << 1., -1.;
  a3 << 1., 0.;
  b1 << 2., 4., 4.;
  b2 << 2., 4., 4.;
  b3 << -2., -10., -7.;
  singleValues(0) = 12.;
  singleValues(1) = 3.;

  // prepare preConditioner to be used to construct a SVD factorization class
  auto prec(std::make_shared<ConstantPreconditioner>(factors));

  // prepare matrix operation to be used in SVD update
  ParallelMatrixOperations matOperation;
  matOperation.initialize(cyclicCommunication);
  auto ptrParMatrixOp = std::make_shared<ParallelMatrixOperations>(matOperation);

  // construct a SVD factorization object
  SVDFactorization svd_1(eps, prec);

  svd_1.initialize(ptrParMatrixOp, m, n);

  svd_1.update(a1, b1);
  svd_1.update(a2, b2);
  svd_1.update(a3, b3);

  // check if the over small single values have been correctly truncated
  BOOST_TEST(svd_1.rank() == 2);
  BOOST_TEST(testing::equals(svd_1.singularValues(), singleValues, 1e-10));
}

BOOST_AUTO_TEST_SUITE_END()

#endif

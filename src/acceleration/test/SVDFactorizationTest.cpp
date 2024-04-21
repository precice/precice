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

BOOST_AUTO_TEST_CASE(testSVDFactorization)
{
  PRECICE_TEST(1_rank);
  int                      m                   = 8;
  double                   eps                 = 1.2e-3;
  const bool               cyclicCommunication = false;
  std::vector<double>      factors(1.0, 1.0);
  Eigen::VectorXd          a(m);
  Eigen::VectorXd          b(m);
  SVDFactorization::Vector singleValues = SVDFactorization::Vector::Zero(2);
  singleValues(0)                       = 0.000578589774323541;
  singleValues(1)                       = 8.13011257500996e-05;

  // prepare preConditioner to be used to construct a SVD factorization class
  auto prec(std::make_shared<impl::ConstantPreconditioner>(factors));

  // prepare matrix operation to be used in SVD update
  ParallelMatrixOperations matOperation;
  matOperation.initialize(cyclicCommunication);
  auto ptrParMatrixOp = std::make_shared<ParallelMatrixOperations>(matOperation);

  // construct a SVD factorization object
  SVDFactorization svd_1(eps, prec);

  svd_1.initialize(ptrParMatrixOp, m, m);
  // update 4 times
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < m; j++) {
      a(j) = j * 0.001 - i * 0.001;
      b(j) = j * 0.001 - i * 0.007;
    }
    svd_1.update(a, b);
  }
  // check if the over small single values have been correctly truncated
  BOOST_TEST(svd_1.rank() == 2);
  BOOST_TEST(testing::equals(svd_1.singularValues(), singleValues, 1e-10));
}

BOOST_AUTO_TEST_SUITE_END()

#endif

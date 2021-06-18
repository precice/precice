#include <Eigen/Core>
#include <algorithm>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "acceleration/MVQNAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice;
using namespace precice::acceleration;

BOOST_AUTO_TEST_SUITE(AccelerationTests)

struct AccelerationSerialTestsFixture {
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;

  //AccelerationSerialTestsFixture() {}
};

BOOST_FIXTURE_TEST_SUITE(AccelerationSerialTests, AccelerationSerialTestsFixture)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(testMVQNPP)
{
  PRECICE_TEST(1_rank);
  //use two vectors and see if underrelaxation works
  double           initialRelaxation        = 0.01;
  int              maxIterationsUsed        = 50;
  int              timestepsReused          = 6;
  int              reusedTimestepsAtRestart = 0;
  int              chunkSize                = 0;
  int              filter                   = Acceleration::QR1FILTER;
  int              restartType              = MVQNAcceleration::NO_RESTART;
  double           singularityLimit         = 1e-10;
  double           svdTruncationEps         = 0.0;
  bool             enforceInitialRelaxation = false;
  bool             alwaysBuildJacobian      = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(factors));
  mesh::PtrMesh           dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));

  MVQNAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                      timestepsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                      restartType, chunkSize, reusedTimestepsAtRestart, svdTruncationEps);

  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;

  //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;

  cplscheme::PtrCouplingData dpcd(new cplscheme::CouplingData(displacements, dummyMesh, false));
  cplscheme::PtrCouplingData fpcd(new cplscheme::CouplingData(forces, dummyMesh, false));

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  forces->values() << 0.1, 0.1, 0.1, 0.1;

  pp.performAcceleration(data);

  BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00000000000000000000));
  BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01000000000000000888));
  BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.02000000000000001776));
  BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.03000000000000002665));
  BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(2), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(3), 0.199000000000000010214));

  data.begin()->second->values() << 10, 10, 10, 10;

  pp.performAcceleration(data);

  BOOST_TEST(testing::equals(data.at(0)->values()(0), -5.63401340929695848558e-01));
  BOOST_TEST(testing::equals(data.at(0)->values()(1), 6.10309919173602111186e-01));
  BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.78402117927690184729e+00));
  BOOST_TEST(testing::equals(data.at(0)->values()(3), 2.95773243938020247157e+00));
  BOOST_TEST(testing::equals(data.at(1)->values()(0), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(1), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(2), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(3), 8.28025852497733250157e-02));
}

BOOST_AUTO_TEST_CASE(testVIQNPP)
{
  PRECICE_TEST(1_rank);
  //use two vectors and see if underrelaxation works

  double           initialRelaxation        = 0.01;
  int              maxIterationsUsed        = 50;
  int              timestepsReused          = 6;
  int              filter                   = acceleration::BaseQNAcceleration::QR1FILTER;
  double           singularityLimit         = 1e-10;
  bool             enforceInitialRelaxation = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  acceleration::impl::PtrPreconditioner prec(new acceleration::impl::ConstantPreconditioner(factors));

  std::map<int, double> scalings;
  scalings.insert(std::make_pair(0, 1.0));
  scalings.insert(std::make_pair(1, 1.0));
  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));

  IQNILSAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                        timestepsReused, filter, singularityLimit, dataIDs, prec);

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;

  //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;

  cplscheme::PtrCouplingData dpcd(new cplscheme::CouplingData(displacements, dummyMesh, false));
  cplscheme::PtrCouplingData fpcd(new cplscheme::CouplingData(forces, dummyMesh, false));
  dpcd->storeIteration();
  fpcd->storeIteration();

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  forces->values() << 0.1, 0.1, 0.1, 0.1;

  pp.performAcceleration(data);

  BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00));
  BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01));
  BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.02));
  BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.03));
  BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->values()(2), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->values()(3), 0.199));

  Eigen::VectorXd newdvalues;
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  data.begin()->second->values() = newdvalues;

  pp.performAcceleration(data);

  BOOST_TEST(testing::equals(data.at(0)->values()(0), -5.63401340929692295845e-01));
  BOOST_TEST(testing::equals(data.at(0)->values()(1), 6.10309919173607440257e-01));
  BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.78402117927690717636e+00));
  BOOST_TEST(testing::equals(data.at(0)->values()(3), 2.95773243938020513610e+00));
  BOOST_TEST(testing::equals(data.at(1)->values()(0), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(1), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(2), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->values()(3), 8.28025852497733944046e-02));
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

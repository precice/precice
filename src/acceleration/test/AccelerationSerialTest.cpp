#include <Eigen/Core>
#include <algorithm>
#include "acceleration/Acceleration.hpp"
#include "acceleration/AitkenAcceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/ConstantRelaxationAcceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "acceleration/IQNIMVJAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/ResidualPreconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "acceleration/test/helper.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice;
using namespace precice::acceleration;

using precice::testing::makeCouplingData;

BOOST_AUTO_TEST_SUITE(AccelerationTests)

struct AccelerationSerialTestsFixture {
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;

  // AccelerationSerialTestsFixture() {}
};

BOOST_FIXTURE_TEST_SUITE(AccelerationSerialTests, AccelerationSerialTestsFixture)

#ifndef PRECICE_NO_MPI

void testIQNIMVJPP(bool exchangeSubsteps)
{
  using DataMap = AccelerationSerialTestsFixture::DataMap;
  // use two vectors and see if underrelaxation works
  double           initialRelaxation          = 0.01;
  int              maxIterationsUsed          = 50;
  int              timeWindowsReused          = 6;
  int              reusedTimeWindowsAtRestart = 0;
  int              chunkSize                  = 0;
  int              filter                     = Acceleration::QR1FILTER;
  int              restartType                = IQNIMVJAcceleration::NO_RESTART;
  double           singularityLimit           = 1e-10;
  double           svdTruncationEps           = 0.0;
  bool             enforceInitialRelaxation   = false;
  bool             alwaysBuildJacobian        = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(factors));
  mesh::PtrMesh           dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));

  IQNIMVJAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                         timeWindowsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                         restartType, chunkSize, reusedTimeWindowsAtRestart, svdTruncationEps, !exchangeSubsteps);

  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  // init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;
  displacements->setSampleAtTime(0, displacements->sample());
  displacements->setSampleAtTime(1, displacements->sample());

  // init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(0, forces->sample());
  forces->setSampleAtTime(1, forces->sample());

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(1, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(1, forces->sample());

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

  // Update the waveform as well
  displacements->setSampleAtTime(1, displacements->sample());
  forces->setSampleAtTime(1, forces->sample());

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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testIQNIMVJPPWithSubsteps)
{
  PRECICE_TEST();
  testIQNIMVJPP(true);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testIQNIMVJPPWithoutSubsteps)
{
  PRECICE_TEST();
  testIQNIMVJPP(false);
}

void testVIQNPP(bool exchangeSubsteps)
{
  using DataMap = AccelerationSerialTestsFixture::DataMap;
  // use two vectors and see if underrelaxation works

  double           initialRelaxation        = 0.01;
  int              maxIterationsUsed        = 50;
  int              timeWindowsReused        = 6;
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
                        timeWindowsReused, filter, singularityLimit, dataIDs, prec, !exchangeSubsteps);

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  // init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;
  displacements->setSampleAtTime(0, displacements->sample());
  displacements->setSampleAtTime(1, displacements->sample());

  // init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(0, forces->sample());
  forces->setSampleAtTime(1, forces->sample());

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  dpcd->storeIteration();
  fpcd->storeIteration();

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(1, displacements->sample());

  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(1, forces->sample());

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

  displacements->setSampleAtTime(1, displacements->sample());
  forces->setSampleAtTime(1, forces->sample());

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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testVIQNPPWithSubsteps)
{
  PRECICE_TEST();
  testVIQNPP(true);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testVIQNPPWithoutSubsteps)
{
  PRECICE_TEST();
  testVIQNPP(false);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithSubsteps)
{
  PRECICE_TEST();
  // use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  mesh::PtrMesh    dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", 3, testing::nextMeshID());

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(0, displacements->sample());

  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(0, forces->sample());

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->setSampleAtTime(1, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(1, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 2);
  BOOST_TEST(data.at(0)->values()(1) == 2);
  BOOST_TEST(data.at(0)->values()(2) == 2.6);
  BOOST_TEST(data.at(0)->values()(3) == 2.8);
  BOOST_TEST(data.at(1)->values()(0) == 0.16);
  BOOST_TEST(data.at(1)->values()(1) == 0.16);
  BOOST_TEST(data.at(1)->values()(2) == 0.16);
  BOOST_TEST(data.at(1)->values()(3) == 0.16);

  displacements->values() << 10, 10, 10, 10;
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 4.6);
  BOOST_TEST(data.at(0)->values()(1) == 5.2);
  BOOST_TEST(data.at(0)->values()(2) == 5.8);
  BOOST_TEST(data.at(0)->values()(3) == 6.4);
  BOOST_TEST(data.at(1)->values()(0) == 0.184);
  BOOST_TEST(data.at(1)->values()(1) == 0.184);
  BOOST_TEST(data.at(1)->values()(2) == 0.184);
  BOOST_TEST(data.at(1)->values()(3) == 0.184);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testAitkenUnderrelaxationWithoutSubsteps)
{
  PRECICE_TEST();

  double              relaxation = 0.4;
  std::vector<int>    dataIDs{0, 1};
  std::vector<double> factors{1, 1};
  mesh::PtrMesh       dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", 3, testing::nextMeshID());

  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(factors));
  AitkenAcceleration      acc(relaxation, dataIDs, prec);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(0.0, displacements->sample());

  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(0.0, forces->sample());

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, false);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, false);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 2);
  BOOST_TEST(data.at(0)->values()(1) == 2);
  BOOST_TEST(data.at(0)->values()(2) == 2.6);
  BOOST_TEST(data.at(0)->values()(3) == 2.8);
  BOOST_TEST(data.at(1)->values()(0) == 0.16);
  BOOST_TEST(data.at(1)->values()(1) == 0.16);
  BOOST_TEST(data.at(1)->values()(2) == 0.16);
  BOOST_TEST(data.at(1)->values()(3) == 0.16);

  data.begin()->second->values() << 10, 10, 10, 10;
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 1.2689851805508461);
  BOOST_TEST(data.at(0)->values()(1) == 2.2390979382674185);
  BOOST_TEST(data.at(0)->values()(2) == 3.2092106959839914);
  BOOST_TEST(data.at(0)->values()(3) == 4.1793234537005644);
  BOOST_TEST(data.at(1)->values()(0) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->values()(1) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->values()(2) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->values()(3) == 0.19880451030866292);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testAitkenUnderrelaxationWithPreconditioner)
{
  PRECICE_TEST();

  double           relaxation = 0.8;
  std::vector<int> dataIDs{0, 1, 2, 3};
  mesh::PtrMesh    dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", 3, testing::nextMeshID());

  impl::PtrPreconditioner prec(new impl::ResidualPreconditioner(-1));
  AitkenAcceleration      acc(relaxation, dataIDs, prec);

  mesh::PtrData data1 = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData data2 = std::make_shared<mesh::Data>("fvalues", -1, 1);
  mesh::PtrData data3 = std::make_shared<mesh::Data>("gvalues", -1, 3);
  mesh::PtrData data4 = std::make_shared<mesh::Data>("hvalues", -1, 1);

  // init data1
  data1->values().resize(2);
  data1->values() << 40, 80;
  data1->setSampleAtTime(0.0, data1->sample());

  // init data2
  data2->values().resize(2);
  data2->values() << 5, 5;
  data2->setSampleAtTime(0.0, data2->sample());

  // init data3
  data3->values().resize(3);
  data3->values() << 1, 2, 3;
  data3->setSampleAtTime(0.0, data3->sample());

  // init data4
  data4->values().resize(4);
  data4->values() << 20, 40, 60, 80;
  data4->setSampleAtTime(0.0, data4->sample());

  cplscheme::PtrCouplingData dpcd = makeCouplingData(data1, dummyMesh, false);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(data2, dummyMesh, false);
  cplscheme::PtrCouplingData gpcd = makeCouplingData(data3, dummyMesh, false);
  cplscheme::PtrCouplingData hpcd = makeCouplingData(data4, dummyMesh, false);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(2, gpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(3, hpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();
  gpcd->storeIteration();
  hpcd->storeIteration();

  acc.initialize(data);

  data1->values() << 1, 7;
  data1->setSampleAtTime(1.0, data1->sample());
  data2->values() << 10, 10;
  data2->setSampleAtTime(1.0, data2->sample());
  data3->values() << 10, 11, 12;
  data3->setSampleAtTime(1.0, data3->sample());
  data4->values() << 40, 60, 80, 100;
  data4->setSampleAtTime(1.0, data4->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 8.8);
  BOOST_TEST(data.at(0)->values()(1) == 21.6);
  BOOST_TEST(data.at(1)->values()(0) == 9);
  BOOST_TEST(data.at(1)->values()(1) == 9);
  BOOST_TEST(data.at(2)->values()(0) == 8.2);
  BOOST_TEST(data.at(2)->values()(1) == 9.2);
  BOOST_TEST(data.at(2)->values()(2) == 10.2);
  BOOST_TEST(data.at(3)->values()(0) == 36);
  BOOST_TEST(data.at(3)->values()(1) == 56);
  BOOST_TEST(data.at(3)->values()(2) == 76);
  BOOST_TEST(data.at(3)->values()(3) == 96);

  data1->values() << 2, 14;
  data1->setSampleAtTime(1.0, data1->sample());
  data2->values() << 8, 8;
  data2->setSampleAtTime(1.0, data2->sample());
  data3->values() << 13, 14, 15;
  data3->setSampleAtTime(1.0, data3->sample());
  data4->values() << 41, 61, 81, 90;
  data4->setSampleAtTime(1.0, data4->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == -17.745640722103754);
  BOOST_TEST(data.at(0)->values()(1) == -20.295060201548626);
  BOOST_TEST(data.at(1)->values()(0) == 9.5588663727976648);
  BOOST_TEST(data.at(1)->values()(1) == 9.5588663727976648);
  BOOST_TEST(data.at(2)->values()(0) == 19.235465491190659);
  BOOST_TEST(data.at(2)->values()(1) == 20.235465491190659);
  BOOST_TEST(data.at(2)->values()(2) == 21.235465491190659);
  BOOST_TEST(data.at(3)->values()(0) == 51.912064609583652);
  BOOST_TEST(data.at(3)->values()(1) == 71.912064609583652);
  BOOST_TEST(data.at(3)->values()(2) == 91.912064609583652);
  BOOST_TEST(data.at(3)->values()(3) == 95.196221242658879);

  data1->values() << 2.1, 14.1;
  data1->setSampleAtTime(1.0, data1->sample());
  data2->values() << 8, 8;
  data2->setSampleAtTime(1.0, data2->sample());
  data3->values() << 13.05, 14.07, 15.1;
  data3->setSampleAtTime(1.0, data3->sample());
  data4->values() << 42, 60, 81.3, 91;
  data4->setSampleAtTime(1.0, data4->sample());

  acc.iterationsConverged(data);

  data1->values() << 3, 16;
  data1->setSampleAtTime(2.0, data1->sample());
  data2->values() << 7, 7;
  data2->setSampleAtTime(2.0, data2->sample());
  data3->values() << 18, 19, 20;
  data3->setSampleAtTime(2.0, data3->sample());
  data4->values() << 50, 70, 90, 110;
  data4->setSampleAtTime(2.0, data4->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 10.4);
  BOOST_TEST(data.at(0)->values()(1) == 28.8);
  BOOST_TEST(data.at(1)->values()(0) == 6.6);
  BOOST_TEST(data.at(1)->values()(1) == 6.6);
  BOOST_TEST(data.at(2)->values()(0) == 14.6);
  BOOST_TEST(data.at(2)->values()(1) == 15.6);
  BOOST_TEST(data.at(2)->values()(2) == 16.6);
  BOOST_TEST(data.at(3)->values()(0) == 44);
  BOOST_TEST(data.at(3)->values()(1) == 64);
  BOOST_TEST(data.at(3)->values()(2) == 84);
  BOOST_TEST(data.at(3)->values()(3) == 104);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradientWithSubsteps)
{
  PRECICE_TEST();
  //use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  const int        dim       = 3;
  mesh::PtrMesh    dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", dim, testing::nextMeshID());

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->requireDataGradient();
  displacements->gradients().resize(dim, 4);
  for (unsigned int r = 0; r < dim; ++r) {
    for (unsigned int c = 0; c < 4; ++c)
      displacements->gradients()(r, c) = r + r * c;
  }
  displacements->setSampleAtTime(0.0, displacements->sample());
  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->requireDataGradient();
  forces->gradients().resize(dim, 4);
  forces->gradients().setConstant(-2);
  forces->setSampleAtTime(0.0, forces->sample());

  bool exchangeSubsteps = true;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->gradients().setConstant(2.5);
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->gradients().setConstant(3);
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  // Test value data
  BOOST_TEST(data.at(0)->values()(0) == 2);
  BOOST_TEST(data.at(0)->values()(1) == 2);
  BOOST_TEST(data.at(0)->values()(2) == 2.6);
  BOOST_TEST(data.at(0)->values()(3) == 2.8);
  BOOST_TEST(data.at(1)->values()(0) == 0.16);
  BOOST_TEST(data.at(1)->values()(1) == 0.16);
  BOOST_TEST(data.at(1)->values()(2) == 0.16);
  BOOST_TEST(data.at(1)->values()(3) == 0.16);

  // Test gradient data
  BOOST_TEST(data.at(0)->gradients()(0, 0) == 1);
  BOOST_TEST(data.at(0)->gradients()(0, 1) == 1);
  BOOST_TEST(data.at(0)->gradients()(0, 2) == 1);
  BOOST_TEST(data.at(0)->gradients()(1, 0) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(1, 1) == 2.2);
  BOOST_TEST(data.at(0)->gradients()(1, 2) == 2.8);
  BOOST_TEST(data.at(1)->gradients()(0, 0) == 0);
  BOOST_TEST(data.at(1)->gradients()(0, 1) == 0);
  BOOST_TEST(data.at(1)->gradients()(0, 2) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 0) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 1) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 2) == 0);

  displacements->values() << 10, 10, 10, 10;
  displacements->gradients().setConstant(4);
  displacements->setSampleAtTime(1, displacements->sample());
  forces->setSampleAtTime(1, forces->sample());

  acc.performAcceleration(data);

  // Check that store iteration works properly
  BOOST_TEST(data.at(0)->values()(0) == 4.6);
  BOOST_TEST(data.at(0)->values()(1) == 5.2);
  BOOST_TEST(data.at(0)->values()(2) == 5.8);
  BOOST_TEST(data.at(0)->values()(3) == 6.4);
  BOOST_TEST(data.at(1)->values()(0) == 0.184);
  BOOST_TEST(data.at(1)->values()(1) == 0.184);
  BOOST_TEST(data.at(1)->values()(2) == 0.184);
  BOOST_TEST(data.at(1)->values()(3) == 0.184);

  BOOST_TEST(data.at(0)->gradients()(0, 0) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(0, 1) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(0, 2) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(1, 0) == 2.2);
  BOOST_TEST(data.at(0)->gradients()(1, 1) == 2.8);
  BOOST_TEST(data.at(0)->gradients()(1, 2) == 3.4);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithoutSubsteps)
{
  PRECICE_TEST();
  //use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  mesh::PtrMesh    dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", 3, testing::nextMeshID());

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(0.0, displacements->sample());

  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(0.0, forces->sample());

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 2);
  BOOST_TEST(data.at(0)->values()(1) == 2);
  BOOST_TEST(data.at(0)->values()(2) == 2.6);
  BOOST_TEST(data.at(0)->values()(3) == 2.8);
  BOOST_TEST(data.at(1)->values()(0) == 0.16);
  BOOST_TEST(data.at(1)->values()(1) == 0.16);
  BOOST_TEST(data.at(1)->values()(2) == 0.16);
  BOOST_TEST(data.at(1)->values()(3) == 0.16);

  displacements->values() << 10, 10, 10, 10;
  displacements->setSampleAtTime(1.0, displacements->sample());
  forces->setSampleAtTime(1.0, forces->sample());

  acc.performAcceleration(data);

  BOOST_TEST(data.at(0)->values()(0) == 4.6);
  BOOST_TEST(data.at(0)->values()(1) == 5.2);
  BOOST_TEST(data.at(0)->values()(2) == 5.8);
  BOOST_TEST(data.at(0)->values()(3) == 6.4);
  BOOST_TEST(data.at(1)->values()(0) == 0.184);
  BOOST_TEST(data.at(1)->values()(1) == 0.184);
  BOOST_TEST(data.at(1)->values()(2) == 0.184);
  BOOST_TEST(data.at(1)->values()(3) == 0.184);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradientWithoutSubsteps)
{
  PRECICE_TEST();
  // use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  const int        dim       = 3;
  mesh::PtrMesh    dummyMesh = std::make_shared<mesh::Mesh>("DummyMesh", dim, testing::nextMeshID());

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->requireDataGradient();
  displacements->gradients().resize(dim, 4);
  for (unsigned int r = 0; r < dim; ++r) {
    for (unsigned int c = 0; c < 4; ++c)
      displacements->gradients()(r, c) = r + r * c;
  }
  displacements->setSampleAtTime(0, displacements->sample());
  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->requireDataGradient();
  forces->gradients().resize(dim, 4);
  forces->gradients().setConstant(-2);
  forces->setSampleAtTime(0, forces->sample());

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->gradients().setConstant(2.5);
  displacements->setSampleAtTime(1, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->gradients().setConstant(3);
  forces->setSampleAtTime(1, forces->sample());

  acc.performAcceleration(data);

  // Test value data
  BOOST_TEST(data.at(0)->values()(0) == 2);
  BOOST_TEST(data.at(0)->values()(1) == 2);
  BOOST_TEST(data.at(0)->values()(2) == 2.6);
  BOOST_TEST(data.at(0)->values()(3) == 2.8);
  BOOST_TEST(data.at(1)->values()(0) == 0.16);
  BOOST_TEST(data.at(1)->values()(1) == 0.16);
  BOOST_TEST(data.at(1)->values()(2) == 0.16);
  BOOST_TEST(data.at(1)->values()(3) == 0.16);

  // Test gradient data
  BOOST_TEST(data.at(0)->gradients()(0, 0) == 1);
  BOOST_TEST(data.at(0)->gradients()(0, 1) == 1);
  BOOST_TEST(data.at(0)->gradients()(0, 2) == 1);
  BOOST_TEST(data.at(0)->gradients()(1, 0) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(1, 1) == 2.2);
  BOOST_TEST(data.at(0)->gradients()(1, 2) == 2.8);
  BOOST_TEST(data.at(1)->gradients()(0, 0) == 0);
  BOOST_TEST(data.at(1)->gradients()(0, 1) == 0);
  BOOST_TEST(data.at(1)->gradients()(0, 2) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 0) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 1) == 0);
  BOOST_TEST(data.at(1)->gradients()(1, 2) == 0);

  data.begin()->second->values() << 10, 10, 10, 10;
  displacements->gradients().setConstant(4);
  displacements->setSampleAtTime(1, displacements->sample());
  forces->setSampleAtTime(1, forces->sample());

  acc.performAcceleration(data);

  // Check that store iteration works properly
  BOOST_TEST(data.at(0)->values()(0) == 4.6);
  BOOST_TEST(data.at(0)->values()(1) == 5.2);
  BOOST_TEST(data.at(0)->values()(2) == 5.8);
  BOOST_TEST(data.at(0)->values()(3) == 6.4);
  BOOST_TEST(data.at(1)->values()(0) == 0.184);
  BOOST_TEST(data.at(1)->values()(1) == 0.184);
  BOOST_TEST(data.at(1)->values()(2) == 0.184);
  BOOST_TEST(data.at(1)->values()(3) == 0.184);

  BOOST_TEST(data.at(0)->gradients()(0, 0) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(0, 1) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(0, 2) == 1.6);
  BOOST_TEST(data.at(0)->gradients()(1, 0) == 2.2);
  BOOST_TEST(data.at(0)->gradients()(1, 1) == 2.8);
  BOOST_TEST(data.at(0)->gradients()(1, 2) == 3.4);
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

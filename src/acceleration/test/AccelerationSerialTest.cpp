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
#include "testing/Meshes.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice;
using namespace precice::acceleration;
using namespace precice::acceleration::impl;

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
  double       initialRelaxation          = 0.01;
  int          maxIterationsUsed          = 50;
  int          timeWindowsReused          = 6;
  int          reusedTimeWindowsAtRestart = 0;
  int          chunkSize                  = 0;
  int          filter                     = Acceleration::QR1FILTER;
  int          restartType                = IQNIMVJAcceleration::NO_RESTART;
  double       singularityLimit           = 1e-10;
  double       svdTruncationEps           = 0.0;
  bool         enforceInitialRelaxation   = false;
  bool         alwaysBuildJacobian        = false;
  const double windowStart                = 0;
  const double windowEnd                  = 1;

  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  PtrPreconditioner prec(new ConstantPreconditioner(factors));
  auto              dummyMesh = testing::makeDummy2DMesh(4);

  IQNIMVJAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                         timeWindowsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                         restartType, chunkSize, reusedTimeWindowsAtRestart, svdTruncationEps, !exchangeSubsteps);

  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  // init displacements & forces
  displacements->emplaceSampleAtTime(windowStart, {1.0, 1.0, 1.0, 1.0});
  displacements->emplaceSampleAtTime(windowEnd, {1.0, 1.0, 1.0, 1.0});
  forces->emplaceSampleAtTime(windowStart, {0.2, 0.2, 0.2, 0.2});
  forces->emplaceSampleAtTime(windowEnd, {0.2, 0.2, 0.2, 0.2});

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  pp.initialize(data);

  displacements->emplaceSampleAtTime(windowEnd, {1.0, 2.0, 3.0, 4.0});
  forces->emplaceSampleAtTime(windowEnd, {0.1, 0.1, 0.1, 0.1});

  pp.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(0), 1.00000000000000000000));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(1), 1.01000000000000000888));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(2), 1.02000000000000001776));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(3), 1.03000000000000002665));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(0), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(1), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(2), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(3), 0.199000000000000010214));

  // Update the waveform as well
  displacements->emplaceSampleAtTime(windowEnd, {10, 10, 10, 10});
  forces->setSampleAtTime(windowEnd, forces->sample());

  pp.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(0), -5.63401340929695848558e-01));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(1), 6.10309919173602111186e-01));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(2), 1.78402117927690184729e+00));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(3), 2.95773243938020247157e+00));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(0), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(1), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(2), 8.28025852497733250157e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(3), 8.28025852497733250157e-02));
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
  const double     windowStart = 0;
  const double     windowEnd   = 1;

  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  PtrPreconditioner prec(new ConstantPreconditioner(factors));

  std::map<int, double> scalings;
  scalings.insert(std::make_pair(0, 1.0));
  scalings.insert(std::make_pair(1, 1.0));
  auto dummyMesh = testing::makeDummy2DMesh(4);

  IQNILSAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                        timeWindowsReused, filter, singularityLimit, dataIDs, prec, !exchangeSubsteps);

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  // init displacements & forces
  displacements->emplaceSampleAtTime(windowStart, {1.0, 1.0, 1.0, 1.0});
  displacements->emplaceSampleAtTime(windowEnd, {1.0, 1.0, 1.0, 1.0});
  forces->emplaceSampleAtTime(windowStart, {0.2, 0.2, 0.2, 0.2});
  forces->emplaceSampleAtTime(windowEnd, {0.2, 0.2, 0.2, 0.2});

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  dpcd->storeIteration();
  fpcd->storeIteration();

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  displacements->emplaceSampleAtTime(windowEnd, {1.0, 2.0, 3.0, 4.0});
  forces->emplaceSampleAtTime(windowEnd, {0.1, 0.1, 0.1, 0.1});

  pp.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(0), 1.00));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(1), 1.01));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(2), 1.02));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(3), 1.03));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(0), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(1), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(2), 0.199));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(3), 0.199));

  Eigen::VectorXd newdvalues;
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);

  displacements->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), newdvalues));
  forces->setSampleAtTime(windowEnd, forces->sample());

  pp.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(0), -5.63401340929692295845e-01));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(1), 6.10309919173607440257e-01));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(2), 1.78402117927690717636e+00));
  BOOST_TEST(testing::equals(data.at(0)->timeStepsStorage().sample(windowEnd)(3), 2.95773243938020513610e+00));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(0), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(1), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(2), 8.28025852497733944046e-02));
  BOOST_TEST(testing::equals(data.at(1)->timeStepsStorage().sample(windowEnd)(3), 8.28025852497733944046e-02));
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
  auto             dummyMesh   = testing::makeDummy3DMesh(4);
  const double     windowStart = 0;
  const double     windowEnd   = 1;

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // init displacements & forces
  displacements->emplaceSampleAtTime(windowStart, {1.0, 2.0, 3.0, 4.0});
  forces->emplaceSampleAtTime(windowStart, {0.2, 0.2, 0.2, 0.2});

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->emplaceSampleAtTime(windowEnd, {3.5, 2.0, 2.0, 1.0});
  forces->emplaceSampleAtTime(windowEnd, {0.1, 0.1, 0.1, 0.1});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 2.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 2.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.16);

  displacements->emplaceSampleAtTime(windowEnd, {10, 10, 10, 10});
  forces->setSampleAtTime(windowEnd, forces->sample());

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 4.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 5.2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 5.8);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 6.4);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.184);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testAitkenUnderrelaxationWithoutSubsteps)
{
  PRECICE_TEST();

  double              relaxation = 0.4;
  std::vector<int>    dataIDs{0, 1};
  std::vector<double> factors{1, 1};
  auto                dummyMesh   = testing::makeDummy3DMesh(4);
  const double        windowStart = 0;
  const double        windowEnd   = 1;

  PtrPreconditioner  prec(new ConstantPreconditioner(factors));
  AitkenAcceleration acc(relaxation, dataIDs, prec);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // //init displacements & forces
  displacements->emplaceSampleAtTime(windowStart, {1.0, 2.0, 3.0, 4.0});
  forces->emplaceSampleAtTime(windowStart, {0.2, 0.2, 0.2, 0.2});

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, false);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, false);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->emplaceSampleAtTime(windowEnd, {3.5, 2.0, 2.0, 1.0});
  forces->emplaceSampleAtTime(windowEnd, {0.1, 0.1, 0.1, 0.1});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 2.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 2.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.16);

  displacements->emplaceSampleAtTime(windowEnd, {10, 10, 10, 10});
  forces->setSampleAtTime(windowEnd, forces->sample());

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 1.2689851805508461);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2.2390979382674185);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 3.2092106959839914);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 4.1793234537005644);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.19880451030866292);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.19880451030866292);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testAitkenUnderrelaxationWithPreconditioner)
{
  PRECICE_TEST();

  double           relaxation = 0.8;
  std::vector<int> dataIDs{0, 1, 2, 3};
  auto             dummyMesh = testing::makeDummy2DMesh(2);

  double       windowStart = 0;
  double       windowEnd   = 1;
  const double dt          = 1;

  PtrPreconditioner  prec(new ResidualPreconditioner(-1));
  AitkenAcceleration acc(relaxation, dataIDs, prec);

  mesh::PtrData data1 = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData data2 = std::make_shared<mesh::Data>("fvalues", -1, 1);
  mesh::PtrData data3 = std::make_shared<mesh::Data>("gvalues", -1, 2);
  mesh::PtrData data4 = std::make_shared<mesh::Data>("hvalues", -1, 2);

  // init data
  data1->emplaceSampleAtTime(windowStart, {40, 80});
  data2->emplaceSampleAtTime(windowStart, {5, 5});
  data3->emplaceSampleAtTime(windowStart, {1, 2, 3, 4});
  data4->emplaceSampleAtTime(windowStart, {20, 40, 60, 80});

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

  data1->emplaceSampleAtTime(windowEnd, {1, 7});
  data2->emplaceSampleAtTime(windowEnd, {10, 10});
  data3->emplaceSampleAtTime(windowEnd, {10, 11, 12, 13});
  data4->emplaceSampleAtTime(windowEnd, {40, 60, 80, 100});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 8.8);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 21.6);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 9);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 9);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(0) == 8.2);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(1) == 9.2);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(2) == 10.2);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(0) == 36);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(1) == 56);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(2) == 76);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(3) == 96);

  data1->emplaceSampleAtTime(windowEnd, {2, 14});
  data2->emplaceSampleAtTime(windowEnd, {8, 8});
  data3->emplaceSampleAtTime(windowEnd, {13, 14, 15, 16});
  data4->emplaceSampleAtTime(windowEnd, {41, 61, 81, 90});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == -17.745640722103754);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == -20.295060201548626);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 9.5588663727976648);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 9.5588663727976648);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(0) == 19.235465491190659);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(1) == 20.235465491190659);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(2) == 21.235465491190659);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(0) == 51.912064609583652);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(1) == 71.912064609583652);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(2) == 91.912064609583652);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(3) == 95.196221242658879);

  data1->emplaceSampleAtTime(windowEnd, {2.1, 14.1});
  data2->emplaceSampleAtTime(windowEnd, {8, 8});
  data3->emplaceSampleAtTime(windowEnd, {13.05, 14.07, 15.1, 16.1});
  data4->emplaceSampleAtTime(windowEnd, {42, 60, 81.3, 91});

  acc.iterationsConverged(data, windowStart);

  // move to next window
  windowStart += dt;
  windowEnd += dt;

  // move to next window
  windowStart += dt;
  windowEnd += dt;

  data1->emplaceSampleAtTime(windowEnd, {3, 16});
  data2->emplaceSampleAtTime(windowEnd, {7, 7});
  data3->emplaceSampleAtTime(windowEnd, {18, 19, 20, 21});
  data4->emplaceSampleAtTime(windowEnd, {50, 70, 90, 110});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 10.4);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 28.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 6.6);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 6.6);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(0) == 14.6);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(1) == 15.6);
  BOOST_TEST(data.at(2)->timeStepsStorage().sample(windowEnd)(2) == 16.6);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(0) == 44);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(1) == 64);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(2) == 84);
  BOOST_TEST(data.at(3)->timeStepsStorage().sample(windowEnd)(3) == 104);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradientWithSubsteps)
{
  PRECICE_TEST();
  // use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  const int        dim         = 3;
  auto             dummyMesh   = testing::makeDummy3DMesh(4);
  const double     windowStart = 0;
  const double     windowEnd   = 1;

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // init displacements
  displacements->requireDataGradient();
  Eigen::MatrixXd displacementGradient(displacements->gradients());
  displacementGradient.resize(dim, 4);
  for (unsigned int r = 0; r < dim; ++r) {
    for (unsigned int c = 0; c < 4; ++c)
      displacementGradient(r, c) = r + r * c;
  }
  displacements->setSampleAtTime(windowStart, time::Sample(displacements->getDimensions(), Eigen::Vector4d{1.0, 2.0, 3.0, 4.0}, displacementGradient));
  // init forces
  forces->requireDataGradient();
  Eigen::MatrixXd forcesGradient(forces->gradients());
  forcesGradient.resize(dim, 4);
  forcesGradient.setConstant(-2);
  forces->setSampleAtTime(windowStart, time::Sample(forces->getDimensions(), Eigen::Vector4d{0.2, 0.2, 0.2, 0.2}, forcesGradient));

  bool exchangeSubsteps = true;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), Eigen::Vector4d{3.5, 2.0, 2.0, 1.0}, Eigen::MatrixXd(displacements->gradients()).setConstant(2.5)));
  forces->setSampleAtTime(windowEnd, time::Sample(forces->getDimensions(), Eigen::Vector4d{0.1, 0.1, 0.1, 0.1}, Eigen::MatrixXd(forces->gradients()).setConstant(3)));

  acc.performAcceleration(data, windowStart, windowEnd);

  // Test value data
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 2.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 2.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.16);

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

  displacements->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), Eigen::Vector4d{10, 10, 10, 10}, Eigen::MatrixXd(displacements->gradients()).setConstant(4)));
  forces->setSampleAtTime(windowEnd, forces->sample());

  acc.performAcceleration(data, windowStart, windowEnd);

  // Check that store iteration works properly
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 4.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 5.2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 5.8);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 6.4);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.184);

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
  // use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  auto             dummyMesh   = testing::makeDummy3DMesh(4);
  const double     windowStart = 0;
  const double     windowEnd   = 1;

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  displacements->emplaceSampleAtTime(windowStart, {1.0, 2.0, 3.0, 4.0});
  forces->emplaceSampleAtTime(windowStart, {0.2, 0.2, 0.2, 0.2});

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->emplaceSampleAtTime(windowEnd, {3.5, 2.0, 2.0, 1.0});
  forces->emplaceSampleAtTime(windowEnd, {0.1, 0.1, 0.1, 0.1});

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 2.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 2.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.16);

  displacements->emplaceSampleAtTime(windowEnd, {10, 10, 10, 10});

  forces->setSampleAtTime(windowEnd, forces->sample());

  acc.performAcceleration(data, windowStart, windowEnd);

  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 4.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 5.2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 5.8);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 6.4);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.184);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradientWithoutSubsteps)
{
  PRECICE_TEST();
  // use two vectors and see if underrelaxation works
  double           relaxation = 0.4;
  std::vector<int> dataIDs{0, 1};
  const int        dim         = 3;
  auto             dummyMesh   = testing::makeDummy3DMesh(4);
  const double     windowStart = 0;
  const double     windowEnd   = 1;

  ConstantRelaxationAcceleration acc(relaxation, dataIDs);

  mesh::PtrData displacements = std::make_shared<mesh::Data>("dvalues", -1, 1);
  mesh::PtrData forces        = std::make_shared<mesh::Data>("fvalues", -1, 1);

  // init displacements
  displacements->requireDataGradient();
  Eigen::MatrixXd displacementGradient(displacements->gradients());
  displacementGradient.resize(dim, 4);
  for (unsigned int r = 0; r < dim; ++r) {
    for (unsigned int c = 0; c < 4; ++c)
      displacementGradient(r, c) = r + r * c;
  }
  displacements->setSampleAtTime(windowStart, time::Sample(displacements->getDimensions(), Eigen::Vector4d{1.0, 2.0, 3.0, 4.0}, displacementGradient));
  // init forces
  forces->requireDataGradient();
  Eigen::MatrixXd forcesGradient(forces->gradients());
  forcesGradient.resize(dim, 4);
  forcesGradient.setConstant(-2);
  forces->setSampleAtTime(windowStart, time::Sample(forces->getDimensions(), Eigen::Vector4d{0.2, 0.2, 0.2, 0.2}, forcesGradient));

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, exchangeSubsteps);
  cplscheme::PtrCouplingData fpcd = makeCouplingData(forces, dummyMesh, exchangeSubsteps);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), Eigen::Vector4d{3.5, 2.0, 2.0, 1.0}, Eigen::MatrixXd(displacements->gradients()).setConstant(2.5)));
  forces->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), Eigen::Vector4d{0.1, 0.1, 0.1, 0.1}, Eigen::MatrixXd(displacements->gradients()).setConstant(3)));

  acc.performAcceleration(data, windowStart, windowEnd);

  // Test value data
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 2.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 2.8);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.16);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.16);

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

  displacements->setSampleAtTime(windowEnd, time::Sample(displacements->getDimensions(), Eigen::Vector4d{10, 10, 10, 10}, Eigen::MatrixXd(displacements->gradients()).setConstant(4)));
  forces->setSampleAtTime(windowEnd, forces->sample());

  acc.performAcceleration(data, windowStart, windowEnd);

  // Check that store iteration works properly
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(0) == 4.6);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(1) == 5.2);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(2) == 5.8);
  BOOST_TEST(data.at(0)->timeStepsStorage().sample(windowEnd)(3) == 6.4);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(0) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(1) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(2) == 0.184);
  BOOST_TEST(data.at(1)->timeStepsStorage().sample(windowEnd)(3) == 0.184);

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

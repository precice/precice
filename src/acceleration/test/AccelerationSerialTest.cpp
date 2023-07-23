#include <Eigen/Core>
#include <algorithm>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/ConstantRelaxationAcceleration.hpp"
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
  double           initialRelaxation          = 0.01;
  int              maxIterationsUsed          = 50;
  int              timeWindowsReused          = 6;
  int              reusedTimeWindowsAtRestart = 0;
  int              chunkSize                  = 0;
  int              filter                     = Acceleration::QR1FILTER;
  int              restartType                = MVQNAcceleration::NO_RESTART;
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

  MVQNAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                      timeWindowsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                      restartType, chunkSize, reusedTimeWindowsAtRestart, svdTruncationEps);

  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_START, displacements->sample());

  //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(time::Storage::WINDOW_START, forces->sample());

  bool exchangeSubsteps = false; // @todo set "true" as soon as acceleration scheme supports subcycling

  cplscheme::PtrCouplingData dpcd(new cplscheme::CouplingData(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER));
  cplscheme::PtrCouplingData fpcd(new cplscheme::CouplingData(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER));

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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
                        timeWindowsReused, filter, singularityLimit, dataIDs, prec);

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_START, displacements->sample());

  //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(time::Storage::WINDOW_START, forces->sample());

  bool exchangeSubsteps = false; // @todo set "true" as soon as acceleration scheme supports subcycling

  cplscheme::PtrCouplingData dpcd(new cplscheme::CouplingData(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER));
  cplscheme::PtrCouplingData fpcd(new cplscheme::CouplingData(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER));
  dpcd->storeIteration();
  fpcd->storeIteration();

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  displacements->values() << 1.0, 2.0, 3.0, 4.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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

BOOST_AUTO_TEST_CASE(testConstantUnderrelaxation)
{
  PRECICE_TEST(1_rank);
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
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());

  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = std::make_shared<cplscheme::CouplingData>(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);
  cplscheme::PtrCouplingData fpcd = std::make_shared<cplscheme::CouplingData>(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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

BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradient)
{
  PRECICE_TEST(1_rank);
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
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->requireDataGradient();
  forces->gradients().resize(dim, 4);
  forces->gradients().setConstant(-2);
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

  bool exchangeSubsteps = false;

  cplscheme::PtrCouplingData dpcd = std::make_shared<cplscheme::CouplingData>(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);
  cplscheme::PtrCouplingData fpcd = std::make_shared<cplscheme::CouplingData>(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->gradients().setConstant(2.5);
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->gradients().setConstant(3);
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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

BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithSubsteps)
{
  PRECICE_TEST(1_rank);
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
  displacements->setSampleAtTime(time::Storage::WINDOW_START, displacements->sample());

  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->setSampleAtTime(time::Storage::WINDOW_START, forces->sample());

  bool exchangeSubsteps = true;

  cplscheme::PtrCouplingData dpcd = std::make_shared<cplscheme::CouplingData>(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);
  cplscheme::PtrCouplingData fpcd = std::make_shared<cplscheme::CouplingData>(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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

BOOST_AUTO_TEST_CASE(testConstantUnderrelaxationWithGradientWithSubsteps)
{
  PRECICE_TEST(1_rank);
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
  displacements->setSampleAtTime(time::Storage::WINDOW_START, displacements->sample());
  // //init forces
  forces->values().resize(4);
  forces->values() << 0.2, 0.2, 0.2, 0.2;
  forces->requireDataGradient();
  forces->gradients().resize(dim, 4);
  forces->gradients().setConstant(-2);
  forces->setSampleAtTime(time::Storage::WINDOW_START, forces->sample());

  bool exchangeSubsteps = true;

  cplscheme::PtrCouplingData dpcd = std::make_shared<cplscheme::CouplingData>(displacements, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);
  cplscheme::PtrCouplingData fpcd = std::make_shared<cplscheme::CouplingData>(forces, dummyMesh, false, exchangeSubsteps, cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(1, fpcd));
  dpcd->storeIteration();
  fpcd->storeIteration();

  acc.initialize(data);

  displacements->values() << 3.5, 2.0, 2.0, 1.0;
  displacements->gradients().setConstant(2.5);
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->values() << 0.1, 0.1, 0.1, 0.1;
  forces->gradients().setConstant(3);
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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
  displacements->setSampleAtTime(time::Storage::WINDOW_END, displacements->sample());
  forces->setSampleAtTime(time::Storage::WINDOW_END, forces->sample());

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

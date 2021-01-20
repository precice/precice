#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "acceleration/MVQNAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

struct ParallelImplicitCouplingSchemeFixture {
  using DataMap = std::map<int, PtrCouplingData>;

  std::string _pathToTests;

  ParallelImplicitCouplingSchemeFixture()
  {
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }
};

BOOST_FIXTURE_TEST_SUITE(ParallelImplicitCouplingSchemeTests, ParallelImplicitCouplingSchemeFixture)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(testParseConfigurationWithRelaxation)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;

  std::string path(_pathToTests + "parallel-implicit-cplscheme-relax-const-config.xml");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, xml::ConfigurationContext{}, path);
  BOOST_CHECK(cplSchemeConfig._accelerationConfig->getAcceleration().get());
}

BOOST_AUTO_TEST_CASE(testMVQNPP)
{
  PRECICE_TEST(1_rank);
  //use two vectors and see if underrelaxation works
  double           initialRelaxation        = 0.01;
  int              maxIterationsUsed        = 50;
  int              timestepsReused          = 6;
  int              reusedTimestepsAtRestart = 0;
  int              chunkSize                = 0;
  int              filter                   = acceleration::Acceleration::QR1FILTER;
  int              restartType              = acceleration::MVQNAcceleration::NO_RESTART;
  double           singularityLimit         = 1e-10;
  double           svdTruncationEps         = 0.0;
  bool             enforceInitialRelaxation = false;
  bool             alwaysBuildJacobian      = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  acceleration::impl::PtrPreconditioner prec(new acceleration::impl::ConstantPreconditioner(factors));
  mesh::PtrMesh                         dummyMesh(new mesh::Mesh("DummyMesh", 3, false, testing::nextMeshID()));

  acceleration::MVQNAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                                    timestepsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                                    restartType, chunkSize, reusedTimestepsAtRestart, svdTruncationEps);

  Eigen::VectorXd dcol1;
  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  utils::append(displacements->values(), 1.0);
  utils::append(displacements->values(), 2.0);
  utils::append(displacements->values(), 3.0);
  utils::append(displacements->values(), 4.0);

  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);

  PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

  //init forces
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);

  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);

  PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

  DataMap data;
  data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  dpcd->oldValues.col(0) = dcol1;
  fpcd->oldValues.col(0) = fcol1;

  pp.performAcceleration(data);

  BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00000000000000000000));
  BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01000000000000000888));
  BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.02000000000000001776));
  BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.03000000000000002665));
  BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(2), 0.199000000000000010214));
  BOOST_TEST(testing::equals(data.at(1)->values()(3), 0.199000000000000010214));

  Eigen::VectorXd newdvalues;
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);

  data.begin()->second->values() = newdvalues;

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
  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, false, testing::nextMeshID()));

  acceleration::IQNILSAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                                      timestepsReused, filter, singularityLimit, dataIDs, prec);

  Eigen::VectorXd dcol1;
  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  //init displacements
  utils::append(displacements->values(), 1.0);
  utils::append(displacements->values(), 2.0);
  utils::append(displacements->values(), 3.0);
  utils::append(displacements->values(), 4.0);

  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);

  PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

  //init forces
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);
  utils::append(forces->values(), 0.1);

  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);

  PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

  DataMap data;
  data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
  data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

  pp.initialize(data);

  dpcd->oldValues.col(0) = dcol1;
  fpcd->oldValues.col(0) = fcol1;

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

BOOST_AUTO_TEST_CASE(testInitializeData)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = true;
  auto m2n                 = context.connectMasters("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();

  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
  const auto    dataID0 = mesh->createData("Data0", 1)->getID();
  const auto    dataID1 = mesh->createData("Data1", 3)->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create a ParallelImplicitCouplingScheme object
  double      maxTime        = 1.0;
  int         maxTimesteps   = 3;
  double      timestepLength = 0.1;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  int         sendDataIndex              = -1;
  int         receiveDataIndex           = -1;
  bool        dataRequiresInitialization = false;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex              = 0;
    receiveDataIndex           = 1;
    dataRequiresInitialization = true;
  } else {
    sendDataIndex              = 1;
    receiveDataIndex           = 0;
    dataRequiresInitialization = true;
  }

  // Create the coupling scheme object
  ParallelCouplingScheme cplScheme(
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
      context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend(mesh->data().at(sendDataIndex), mesh, dataRequiresInitialization);
  cplScheme.addDataToReceive(mesh->data().at(receiveDataIndex), mesh, dataRequiresInitialization);

  // Add convergence measures
  int                                    minIterations = 3;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure2(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(mesh->data().at(1), false, false, minIterationConvMeasure1, true);
  cplScheme.addConvergenceMeasure(mesh->data().at(0), false, false, minIterationConvMeasure2, true);

  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  cplScheme.initialize(0.0, 0);

  if (context.isNamed(nameParticipant0)) {
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    mesh->data(dataID0)->values() = Eigen::VectorXd::Constant(1, 4.0);
    cplScheme.markActionFulfilled(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    auto &values = mesh->data(dataID1)->values();
    BOOST_TEST(testing::equals(values, Eigen::Vector3d(1.0, 2.0, 3.0)), values);

    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(writeIterationCheckpoint)) {
        cplScheme.markActionFulfilled(writeIterationCheckpoint);
      }
      if (cplScheme.isActionRequired(readIterationCheckpoint)) {
        cplScheme.markActionFulfilled(readIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
    }
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    auto &values = mesh->data(dataID0)->values();
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    mesh->data(dataID1)->values() = v;
    cplScheme.markActionFulfilled(constants::actionWriteInitialData());
    BOOST_TEST(testing::equals(values(0), 0.0), values);
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(values(0), 4.0), values);

    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(writeIterationCheckpoint)) {
        cplScheme.markActionFulfilled(writeIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
      if (cplScheme.isActionRequired(readIterationCheckpoint)) {
        cplScheme.markActionFulfilled(readIterationCheckpoint);
      }
    }
  }
  cplScheme.finalize();
}
#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "acceleration/MVQNAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "acceleration/impl/ResidualSumPreconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice;
using namespace precice::cplscheme;
using namespace precice::acceleration;
using namespace precice::acceleration::impl;

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using DataMap = std::map<int, PtrCouplingData>;

BOOST_AUTO_TEST_SUITE(AccelerationMasterSlaveTests)

/// Test that runs on 4 processors.
BOOST_AUTO_TEST_CASE(testVIQNILSpp)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  double           initialRelaxation        = 0.01;
  int              maxIterationsUsed        = 50;
  int              timestepsReused          = 6;
  int              filter                   = BaseQNAcceleration::QR1FILTER;
  double           singularityLimit         = 1e-10;
  bool             enforceInitialRelaxation = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  PtrPreconditioner prec(new ConstantPreconditioner(factors));
  std::vector<int>  vertexOffsets{4, 8, 8, 10};

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, false, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  IQNILSAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                        timestepsReused, filter, singularityLimit, dataIDs, prec);

  Eigen::VectorXd dcol1;
  Eigen::VectorXd fcol1;

  DataMap       data;
  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  if (context.isMaster()) { //Master
    /**
     * processor with 4 vertices
     */

    //init displacements
    Eigen::VectorXd insert(4);
    insert << 1.0, 2.0, 3.0, 4.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(1)) { //Slave1

    /**
     * processor with 4 vertices
     */

    //init displacements
    Eigen::VectorXd insert(4);
    insert << 5.0, 6.0, 7.0, 8.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(2)) { //Slave2

    /**
     * processor with no vertices
     */

    //init displacements
    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(3)) { //Slave3

    /**
     * processor with 2 vertices
     */

    //init displacements
    Eigen::VectorXd insert(2);
    insert << 1.0, 2.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;
  }

  pp.performAcceleration(data);

  Eigen::VectorXd newdvalues;
  if (context.isMaster()) { //Master

    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.02), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.03), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 0.199), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 0.199), data.at(1)->values()(3));

    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);

  } else if (context.isRank(1)) { //Slave1

    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.04), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.05), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.06), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.07), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 0.199), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 0.199), data.at(1)->values()(3));

    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);

  } else if (context.isRank(2)) { //Slave2
    // empty proc
  } else if (context.isRank(3)) { //Slave3

    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 0.199), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 0.199), data.at(1)->values()(1));

    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
  }

  data.begin()->second->values() = newdvalues;

  pp.performAcceleration(data);

  if (context.isMaster()) { //Master
    BOOST_TEST(testing::equals(data.at(0)->values()(0), -1.51483105223442748866e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), -2.35405379763935940218e-01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.04402029270655560822e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 2.32344596517704804484e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254212854123e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254212854123e-02), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 7.23368584254212854123e-02), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 7.23368584254212854123e-02), data.at(1)->values()(3));
  } else if (context.isRank(1)) { //Slave1
    BOOST_TEST(testing::equals(data.at(0)->values()(0), 3.60287163764754048145e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 4.88229731011803202989e+00), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 6.16172298258852357833e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 7.44114865505901601495e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254212854123e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254212854123e-02), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 7.23368584254212854123e-02), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 7.23368584254212854123e-02), data.at(1)->values()(3));
  } else if (context.isRank(2)) { //Slave2
    // empty proc
  } else if (context.isRank(3)) { //Slave3
    BOOST_TEST(testing::equals(data.at(0)->values()(0), -1.51483105223442748866e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), -2.35405379763935940218e-01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254212854123e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254212854123e-02), data.at(1)->values()(1));
  }
}

/// Test that runs on 4 processors.
BOOST_AUTO_TEST_CASE(testVIQNIMVJpp)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  double initialRelaxation        = 0.01;
  int    maxIterationsUsed        = 50;
  int    timestepsReused          = 6;
  int    filter                   = BaseQNAcceleration::QR1FILTER;
  int    restartType              = MVQNAcceleration::NO_RESTART;
  int    chunkSize                = 0;
  int    reusedTimeStepsAtRestart = 0;
  double singularityLimit         = 1e-10;
  double svdTruncationEps         = 0.0;
  bool   enforceInitialRelaxation = false;
  bool   alwaysBuildJacobian      = false;

  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2, 1.0);
  PtrPreconditioner prec(new ConstantPreconditioner(factors));
  std::vector<int>  vertexOffsets{4, 8, 8, 10};

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, false, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  MVQNAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                      timestepsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                      restartType, chunkSize, reusedTimeStepsAtRestart, svdTruncationEps);

  Eigen::VectorXd dcol1;
  Eigen::VectorXd fcol1;

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 1));

  DataMap data;

  if (context.isMaster()) { //Master

    /**
     * processor with 4 vertices
     */

    //init displacements
    Eigen::VectorXd insert(4);
    insert << 1.0, 2.0, 3.0, 4.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(1)) { //Slave1

    /**
     * processor with 4 vertices
     */

    //init displacements
    Eigen::VectorXd insert(4);
    insert << 5.0, 6.0, 7.0, 8.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(2)) { //Slave2

    /**
     * processor with no vertices
     */

    //init displacements
    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;

  } else if (context.isRank(3)) { //Slave3

    /**
     * processor with 2 vertices
     */

    //init displacements
    Eigen::VectorXd insert(2);
    insert << 1.0, 2.0;
    utils::append(displacements->values(), insert);
    insert << 1.0, 1.0;
    utils::append(dcol1, insert);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    insert << 0.1, 0.1;
    utils::append(forces->values(), insert);
    insert << 0.2, 0.2;
    utils::append(fcol1, insert);

    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(1, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = dcol1;
    fpcd->oldValues.col(0) = fcol1;
  }

  pp.performAcceleration(data);

  Eigen::VectorXd newdvalues;
  if (context.isMaster()) { //Master
    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00000000000000000000e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01000000000000000888e+00), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.02000000000000001776e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.03000000000000002665e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 1.99000000000000010214e-01), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 1.99000000000000010214e-01), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 1.99000000000000010214e-01), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 1.99000000000000010214e-01), data.at(1)->values()(3));
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
  } else if (context.isRank(1)) { //Slave1
    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.04000000000000003553e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.05000000000000004441e+00), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.06000000000000005329e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 1.07000000000000006217e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 1.99000000000000010214e-01), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 1.99000000000000010214e-01), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 1.99000000000000010214e-01), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 1.99000000000000010214e-01), data.at(1)->values()(3));
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
  } else if (context.isRank(2)) { //Slave2
    // empty proc
  } else if (context.isRank(3)) { //Slave3
    BOOST_TEST(testing::equals(data.at(0)->values()(0), 1.00000000000000000000e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 1.01000000000000000888e+00), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 1.99000000000000010214e-01), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 1.99000000000000010214e-01), data.at(1)->values()(1));
    utils::append(newdvalues, 10.0);
    utils::append(newdvalues, 10.0);
  }

  data.begin()->second->values() = newdvalues;
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master
    BOOST_TEST(testing::equals(data.at(0)->values()(0), -1.51483105223442748866e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), -2.35405379763935940218e-01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 1.04402029270655738458e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 2.32344596517704893301e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254213131679e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254213131679e-02), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 7.23368584254213131679e-02), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 7.23368584254213131679e-02), data.at(1)->values()(3));
  } else if (context.isRank(1)) { //Slave1
    BOOST_TEST(testing::equals(data.at(0)->values()(0), 3.60287163764754048145e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), 4.88229731011803202989e+00), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(0)->values()(2), 6.16172298258852446651e+00), data.at(0)->values()(2));
    BOOST_TEST(testing::equals(data.at(0)->values()(3), 7.44114865505901601495e+00), data.at(0)->values()(3));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254213131679e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254213131679e-02), data.at(1)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(2), 7.23368584254213131679e-02), data.at(1)->values()(2));
    BOOST_TEST(testing::equals(data.at(1)->values()(3), 7.23368584254213131679e-02), data.at(1)->values()(3));
  } else if (context.isRank(2)) { //Slave2
    // empty proc
  } else if (context.isRank(3)) { //Slave3
    BOOST_TEST(testing::equals(data.at(0)->values()(0), -1.51483105223442748866e+00), data.at(0)->values()(0));
    BOOST_TEST(testing::equals(data.at(0)->values()(1), -2.35405379763935940218e-01), data.at(0)->values()(1));
    BOOST_TEST(testing::equals(data.at(1)->values()(0), 7.23368584254213131679e-02), data.at(1)->values()(0));
    BOOST_TEST(testing::equals(data.at(1)->values()(1), 7.23368584254213131679e-02), data.at(1)->values()(1));
  }
}

/// Test that runs on 4 processors.
BOOST_AUTO_TEST_CASE(testIMVJ_effUpdate_pp)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  // config:
  double initialRelaxation        = 0.1;
  int    maxIterationsUsed        = 30;
  int    timestepsReused          = 0;
  int    filter                   = BaseQNAcceleration::QR2FILTER;
  int    restartType              = MVQNAcceleration::NO_RESTART;
  int    chunkSize                = 0;
  int    reusedTimeStepsAtRestart = 0;
  double singularityLimit         = 1e-2;
  double svdTruncationEps         = 0.0;
  bool   enforceInitialRelaxation = false;
  bool   alwaysBuildJacobian      = false;

  std::vector<int> dataIDs;
  dataIDs.push_back(4);
  dataIDs.push_back(5);
  PtrPreconditioner _preconditioner = PtrPreconditioner(new ResidualSumPreconditioner(-1));
  std::vector<int>  vertexOffsets{0, 11, 22, 22};

  mesh::PtrMesh dummyMesh(new mesh::Mesh("dummyMesh", 2, false, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  MVQNAcceleration pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                      timestepsReused, filter, singularityLimit, dataIDs, _preconditioner, alwaysBuildJacobian,
                      restartType, chunkSize, reusedTimeStepsAtRestart, svdTruncationEps);

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 2));
  mesh::PtrData forces(new mesh::Data("fvalues", -1, 2));

  Eigen::VectorXd doldValues;
  Eigen::VectorXd foldValues;

  Eigen::VectorXd dref;
  Eigen::VectorXd fref;
  double          drefNorm = 0., frefNorm = 0.;

  DataMap data;

  if (context.isMaster()) { //Master
    /**
     * processor with no vertices
     */

    displacements->values().resize(0);
    doldValues.resize(0);
    forces->values().resize(0);
    foldValues.resize(0);

    //init displacements
    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));

    //init forces
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(4, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(5, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = doldValues;
    fpcd->oldValues.col(0) = foldValues;

  } else if (context.isRank(1)) { //Slave1
    /**
     * processor with 4 vertices
     */

    //init displacements
    displacements->values().resize(22);
    doldValues.resize(22);
    forces->values().resize(22);
    foldValues.resize(22);
    displacements->values() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    doldValues << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    forces->values() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    foldValues << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(4, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(5, fpcd));

    pp.initialize(data);

    forces->values() << -0.01765744149520443, -0.000534499502588083, 0.05397520666020422, 0.0005546984205735067, 0.05213823386543703, 0.0007618478879228568, -0.01944857239806249, -0.0009206665792022876, -0.02459872346309381, -0.001296931976456198, 0.04688718434761113, 0.001346643628716769, -0.01063536095060684, -0.01905148710330257, 0.02514593936525903, -0.01643393169986981, -0.02189723835016068, -0.000912218689367709, 0.04985117008772211, 0.0009615805506705544, 0.05534647415570375, 0.0004068469082890895;

    dpcd->oldValues.col(0) = doldValues;
    fpcd->oldValues.col(0) = foldValues;

  } else if (context.isRank(2)) { //Slave2
    /**
     * processor with 4 vertices
     */

    //init displacements
    displacements->values().resize(22);
    doldValues.resize(22);
    forces->values().resize(22);
    foldValues.resize(22);
    displacements->values() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    doldValues << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    forces->values() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    foldValues << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(4, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(5, fpcd));

    pp.initialize(data);

    forces->values() << -0.01465589151503364, -0.0002670111835650672, 0.05711438689366102, 0.0002383730129847531, -0.01536098575916998, -0.000285287812066552, 0.05638274807579218, 0.000283961973555227, -0.006856432131857973, -0.006815594391460808, 0.02901925611525407, -0.02907380915674757, 0.05800715138289463, 9.667376010126116e-05, -0.01376443700165205, -9.547563271960956e-05, 0.05768190311116184, 0.0001311583226994801, -0.01408147387131287, -0.0001216961377915992, -0.0163823504288376, -0.0003874626690545313;

    dpcd->oldValues.col(0) = doldValues;
    fpcd->oldValues.col(0) = foldValues;
  } else if (context.isRank(3)) { //Slave3
    /**
     * processor with no vertices
     */

    displacements->values().resize(0);
    doldValues.resize(0);
    forces->values().resize(0);
    foldValues.resize(0);

    PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));
    PtrCouplingData fpcd(new CouplingData(forces, dummyMesh, false));

    data.insert(std::pair<int, PtrCouplingData>(4, dpcd));
    data.insert(std::pair<int, PtrCouplingData>(5, fpcd));

    pp.initialize(data);

    dpcd->oldValues.col(0) = doldValues;
    fpcd->oldValues.col(0) = foldValues;
  }

  // underrelaxation, first iteration
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master

  } else if (context.isRank(1)) { //Slave1

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    fref << -0.001765744149520443, -5.344995025880831e-05, 0.005397520666020422, 5.546984205735068e-05, 0.005213823386543704, 7.618478879228569e-05, -0.001944857239806249, -9.206665792022876e-05, -0.002459872346309381, -0.0001296931976456198, 0.004688718434761114, 0.000134664362871677, -0.001063536095060684, -0.001905148710330257, 0.002514593936525903, -0.001643393169986981, -0.002189723835016068, -9.12218689367709e-05, 0.004985117008772211, 9.615805506705544e-05, 0.005534647415570375, 4.068469082890896e-05;
    frefNorm = 0.01286041129960619;
    drefNorm = 0.0;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.790053057185293e-06, -2.44566429072041e-08, 1.889281703254964e-06, -1.972492834475447e-07, 1.681634609242917e-06, -2.373356532433882e-07, 1.585003447958184e-06, -5.301850772916681e-08, 1.274187257620066e-06, -2.137488936999111e-07, 1.362955262700412e-06, -2.762153471191986e-07, 1.249747540920782e-06, -3.196338173465977e-07, 1.333501893726392e-06, -3.161541101487353e-07, 1.394538527892028e-06, -1.166536323805688e-07, 1.488382850875808e-06, -2.605379508545059e-07, 2.056077021837937e-06, -1.341692715765341e-07;
    forces->values() << -0.01765744187144705, -0.000534499502451157, 0.05397520721069472, 0.0005546984181272257, 0.05213823442800309, 0.000761847881060882, -0.01944857277019029, -0.0009206665773591249, -0.02459872381892192, -0.001296931982922439, 0.04688718490326162, 0.001346643636856003, -0.01063536111298416, -0.01905148734069537, 0.02514593966043068, -0.01643393192020026, -0.02189723869781963, -0.0009122186870252733, 0.04985117065739809, 0.0009615805515004192, 0.05534647470527156, 0.0004068469091761907;
    foldValues << -0.001765744149520443, -5.344995025880831e-05, 0.005397520666020422, 5.546984205735068e-05, 0.005213823386543704, 7.618478879228569e-05, -0.001944857239806249, -9.206665792022876e-05, -0.002459872346309381, -0.0001296931976456198, 0.004688718434761114, 0.000134664362871677, -0.001063536095060684, -0.001905148710330257, 0.002514593936525903, -0.001643393169986981, -0.002189723835016068, -9.12218689367709e-05, 0.004985117008772211, 9.615805506705544e-05, 0.005534647415570375, 4.068469082890896e-05;

  } else if (context.isRank(2)) { //Slave2

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    fref << -0.001465589151503364, -2.670111835650672e-05, 0.005711438689366103, 2.383730129847531e-05, -0.001536098575916998, -2.85287812066552e-05, 0.005638274807579218, 2.83961973555227e-05, -0.0006856432131857974, -0.0006815594391460808, 0.002901925611525407, -0.002907380915674757, 0.005800715138289463, 9.667376010126117e-06, -0.001376443700165206, -9.547563271960956e-06, 0.005768190311116184, 1.311583226994801e-05, -0.001408147387131287, -1.216961377915992e-05, -0.00163823504288376, -3.874626690545313e-05;
    frefNorm = 0.01265754337961098;
    drefNorm = 0.0;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.848184969639987e-06, -1.983566187932991e-07, 1.952383060128974e-06, 1.050101286643166e-07, 2.020975712018586e-06, -9.297459906882382e-08, 2.123910878481957e-06, -3.349554682884977e-08, 0, 0, 0, 0, 7.715047421278781e-07, 2.958323850532032e-07, 6.5137785527863e-07, -3.40165313149562e-07, 1.498023570500414e-06, 2.492038233690158e-07, 1.395223018993416e-06, -3.150663149441921e-07, 1.954718171910318e-06, -3.415637300374603e-08;
    forces->values() << -0.0146558918972568, -0.000267011181975166, 0.05711438744699839, 0.0002383730136872111, -0.0153609861368436, -0.0002852878106683293, 0.05638274862725741, 0.0002839619744993407, -0.00685643232676097, -0.006815594586569211, 0.02901925639144463, -0.02907380943293575, 0.05800715193585099, 9.667375963025685e-05, -0.01376443739049903, -9.547563172575954e-05, 0.05768190366530584, 0.0001311583223016465, -0.01408147425699792, -0.0001216961368213471, -0.01638235080508845, -0.0003874626694560972;
    foldValues << -0.001465589151503364, -2.670111835650672e-05, 0.005711438689366103, 2.383730129847531e-05, -0.001536098575916998, -2.85287812066552e-05, 0.005638274807579218, 2.83961973555227e-05, -0.0006856432131857974, -0.0006815594391460808, 0.002901925611525407, -0.002907380915674757, 0.005800715138289463, 9.667376010126117e-06, -0.001376443700165206, -9.547563271960956e-06, 0.005768190311116184, 1.311583226994801e-05, -0.001408147387131287, -1.216961377915992e-05, -0.00163823504288376, -3.874626690545313e-05;
  } else if (context.isRank(3)) { //Slave3
    // Dummy Slave to be able to reuse the 4 proc Master Slave Fixture
  }

  data.at(4)->oldValues.col(0) = doldValues;
  data.at(5)->oldValues.col(0) = foldValues;

  // QN- Update, 2. iteration
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master

  } else if (context.isRank(1)) { //Slave1

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 2.182991240484406e-07, -2.982517027848927e-09, 2.30400176824822e-07, -2.405478743936717e-08, 2.050773638768581e-07, -2.89433684663868e-08, 1.932930774951963e-07, -6.465670807451687e-09, 1.553886691210811e-07, -2.606693476135542e-08, 1.662140341429634e-07, -3.368479391313121e-08, 1.524082162646531e-07, -3.897972859683727e-08, 1.62622160359393e-07, -3.85553741174048e-08, 1.700656513328811e-07, -1.422605082208521e-08, 1.815100794307223e-07, -3.17729165761968e-08, 2.507410666078866e-07, -1.636210409619326e-08;
    fref << -0.01765744154108767, -0.0005344995025713848, 0.0539752067273372, 0.0005546984202751798, 0.05213823393404264, 0.0007618478870860308, -0.01944857244344392, -0.0009206665789775117, -0.02459872350648748, -0.001296931977244764, 0.04688718441537337, 0.001346643629709359, -0.01063536097040895, -0.01905148713225291, 0.02514593940125557, -0.01643393172673937, -0.0218972383925581, -0.0009122186890820461, 0.04985117015719479, 0.0009615805507717574, 0.05534647422272421, 0.0004068469083972726;
    drefNorm = 6.435143632392166e-07;
    frefNorm = 0.1286041131776192;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.790034504773721e-05, -2.446591076368466e-07, 1.889267115021718e-05, -1.972643201602028e-06, 1.681613350812527e-05, -2.373460013995369e-06, 1.584978895355817e-05, -5.302446869164338e-07, 1.274157692078479e-05, -2.137546278211264e-06, 1.362926508984742e-05, -2.762211725309514e-06, 1.249719424608544e-05, -3.19640295598053e-06, 1.333474052315949e-05, -3.16159193819195e-06, 1.394510078525391e-05, -1.166587691625877e-06, 1.488356439901566e-05, -2.605456452904905e-06, 2.056070000286195e-05, -1.341920935569228e-06;
    doldValues << 2.182991240484406e-07, -2.982517027848927e-09, 2.30400176824822e-07, -2.405478743936717e-08, 2.050773638768581e-07, -2.89433684663868e-08, 1.932930774951963e-07, -6.465670807451687e-09, 1.553886691210811e-07, -2.606693476135542e-08, 1.662140341429634e-07, -3.368479391313121e-08, 1.524082162646531e-07, -3.897972859683727e-08, 1.62622160359393e-07, -3.85553741174048e-08, 1.700656513328811e-07, -1.422605082208521e-08, 1.815100794307223e-07, -3.17729165761968e-08, 2.507410666078866e-07, -1.636210409619326e-08;
    forces->values() << -0.01845221261910751, -0.000473660279052688, 0.05115217408647257, 0.0004997712466226923, 0.04953299847638116, 0.0006834967349236343, -0.02003343430934222, -0.000812620519068711, -0.02461452082338934, -0.00115320030035888, 0.04486850176987132, 0.00121794049154945, -0.01080003301957858, -0.01839753254507924, 0.02398381340216066, -0.01606849579606463, -0.02220297183992202, -0.0008082032560853196, 0.04750952330271016, 0.0008672013640382038, 0.05236940113731539, 0.0003710183715860552;
    foldValues << -0.01765744154108767, -0.0005344995025713848, 0.0539752067273372, 0.0005546984202751798, 0.05213823393404264, 0.0007618478870860308, -0.01944857244344392, -0.0009206665789775117, -0.02459872350648748, -0.001296931977244764, 0.04688718441537337, 0.001346643629709359, -0.01063536097040895, -0.01905148713225291, 0.02514593940125557, -0.01643393172673937, -0.0218972383925581, -0.0009122186890820461, 0.04985117015719479, 0.0009615805507717574, 0.05534647422272421, 0.0004068469083972726;

  } else if (context.isRank(2)) { //Slave2

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 2.253883807144274e-07, -2.418982831708627e-08, 2.380954632167933e-07, 1.280611153486133e-08, 2.464604196428374e-07, -1.133836421999323e-08, 2.590134870407767e-07, -4.084822236363769e-09, 0, 0, 0, 0, 9.408593154806163e-08, 3.607711529166728e-08, 7.943631316463149e-08, -4.148356921273439e-08, 1.826857767882953e-07, 3.039070609254422e-08, 1.701491258462493e-07, -3.842271618341636e-08, 2.38380232908047e-07, -4.165410783473636e-09;
    fref << -0.01465589156164622, -0.0002670111833711768, 0.05711438696114118, 0.0002383730130704187, -0.01536098580522773, -0.0002852878118960371, 0.05638274814304403, 0.0002839619736703628, -0.006856432155626628, -0.006815594415254513, 0.02901925614893584, -0.02907380919042905, 0.05800715145032832, 9.667376004382162e-05, -0.01376443704907241, -9.547563259840837e-05, 0.05768190317874038, 0.0001311583226509638, -0.01408147391834763, -0.0001216961376732758, -0.01638235047472185, -0.0003874626691035027;
    drefNorm = 6.131610103923933e-07;
    frefNorm = 0.1265754339635572;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.848182952307335e-05, -1.983938722952872e-06, 1.952389995095743e-05, 1.049689886611777e-06, 2.020972044646931e-05, -9.30012125294331e-07, 2.123911759834233e-05, -3.352823479948144e-07, 0, 0, 0, 0, 7.715124780435689e-06, 2.958056858428718e-06, 6.513639301665504e-06, -3.401886529062288e-06, 1.498034283416962e-05, 2.491634858078641e-06, 1.39521486945152e-05, -3.151050708450101e-06, 1.954707223943552e-05, -3.417246252999375e-07;
    doldValues << 2.253883807144274e-07, -2.418982831708627e-08, 2.380954632167933e-07, 1.280611153486133e-08, 2.464604196428374e-07, -1.133836421999323e-08, 2.590134870407767e-07, -4.084822236363769e-09, 0, 0, 0, 0, 9.408593154806163e-08, 3.607711529166728e-08, 7.943631316463149e-08, -4.148356921273439e-08, 1.826857767882953e-07, 3.039070609254422e-08, 1.701491258462493e-07, -3.842271618341636e-08, 2.38380232908047e-07, -4.165410783473636e-09;
    forces->values() << -0.01568208277628194, -0.0002595395446636614, 0.0540328986967421, 0.0002362571305830931, -0.01637736854863682, -0.0002699645831085989, 0.05331751790879287, 0.0002707054191427001, -0.007277539612331946, -0.007235194100552225, 0.02757151633202504, -0.02762772092892902, 0.05505877464319012, 0.0001052840945529276, -0.01465499974491537, -0.0001017767294585529, 0.05464614037258596, 0.0001424559420056945, -0.01506072500921042, -0.0001315030046882618, -0.0173164149989076, -0.0003474184175392483;
    foldValues << -0.01465589156164622, -0.0002670111833711768, 0.05711438696114118, 0.0002383730130704187, -0.01536098580522773, -0.0002852878118960371, 0.05638274814304403, 0.0002839619736703628, -0.006856432155626628, -0.006815594415254513, 0.02901925614893584, -0.02907380919042905, 0.05800715145032832, 9.667376004382162e-05, -0.01376443704907241, -9.547563259840837e-05, 0.05768190317874038, 0.0001311583226509638, -0.01408147391834763, -0.0001216961376732758, -0.01638235047472185, -0.0003874626691035027;
  } else if (context.isRank(3)) { //Slave3
    // Dummy Slave to be able to reuse the 4 proc Master Slave Fixture
  }

  data.at(4)->oldValues.col(0) = doldValues;
  data.at(5)->oldValues.col(0) = foldValues;

  // QN- Update, 3. iteration
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master

  } else if (context.isRank(1)) { //Slave1

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.717512702194643e-05, -2.347247264125579e-07, 1.812723818498285e-05, -1.892683301884884e-06, 1.613485037277307e-05, -2.277271314691775e-06, 1.520766651015439e-05, -5.087470614753933e-07, 1.222540039026145e-05, -2.050926769187082e-06, 1.30771205492404e-05, -2.650282676556395e-06, 1.199091589196657e-05, -3.066880426987482e-06, 1.279452699120483e-05, -3.033483087437419e-06, 1.338015894263769e-05, -1.11930952734218e-06, 1.428059443561221e-05, -2.499874235762116e-06, 1.972766650768916e-05, -1.287497600071812e-06;
    fref << -0.01823457911070198, -0.0004903200524407003, 0.05192521422826112, 0.0005148121016499705, 0.05024639863781041, 0.0007049518329474636, -0.01987328081455485, -0.0008422070677716292, -0.02461019581850909, -0.001192558759870676, 0.04542128452828625, 0.001253183648371612, -0.01075494078801318, -0.01857660729355974, 0.02430204227390975, -0.01616856463954442, -0.02211925279965755, -0.0008366860858679643, 0.04815074428053346, 0.0008930454810143837, 0.05318462260012811, 0.0003808294014608455;
    drefNorm = 5.062970166651817e-05;
    frefNorm = 0.1247902554601672;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.659080663925766e-05, -2.839283676791931e-07, 1.756292801739508e-05, -1.881726812992964e-06, 1.564798101471437e-05, -2.265931706775091e-06, 1.470124517392331e-05, -5.705378156142171e-07, 1.186047603634431e-05, -2.115667271562722e-06, 1.273027556448604e-05, -2.674541973319838e-06, 1.165645170777486e-05, -3.135385366949176e-06, 1.247728214631633e-05, -3.082564251671268e-06, 1.295443089215965e-05, -1.185450561958201e-06, 1.387356342346108e-05, -2.500933334689963e-06, 1.911143938064833e-05, -1.289577439500651e-06;
    doldValues << 1.717512702194643e-05, -2.347247264125579e-07, 1.812723818498285e-05, -1.892683301884884e-06, 1.613485037277307e-05, -2.277271314691775e-06, 1.520766651015439e-05, -5.087470614753933e-07, 1.222540039026145e-05, -2.050926769187082e-06, 1.30771205492404e-05, -2.650282676556395e-06, 1.199091589196657e-05, -3.066880426987482e-06, 1.279452699120483e-05, -3.033483087437419e-06, 1.338015894263769e-05, -1.11930952734218e-06, 1.428059443561221e-05, -2.499874235762116e-06, 1.972766650768916e-05, -1.287497600071812e-06;
    forces->values() << -0.08018882094302014, 0.004252226533647444, -0.1681454949889014, -0.00376664315520894, -0.1528473285523723, -0.005402482563102148, -0.06546415562528343, 0.007580544006107058, -0.02584070375141325, 0.01001310706105286, -0.1119519940946071, -0.008779400784032571, -0.02359260557541131, 0.03241153400001811, -0.06629773012631451, 0.01232457880930212, -0.04595115217779744, 0.007271575191182346, -0.1343966169086029, -0.006463593204557456, -0.1788931815052193, -0.00241194798896828;
    foldValues << -0.01823457911070198, -0.0004903200524407003, 0.05192521422826112, 0.0005148121016499705, 0.05024639863781041, 0.0007049518329474636, -0.01987328081455485, -0.0008422070677716292, -0.02461019581850909, -0.001192558759870676, 0.04542128452828625, 0.001253183648371612, -0.01075494078801318, -0.01857660729355974, 0.02430204227390975, -0.01616856463954442, -0.02211925279965755, -0.0008366860858679643, 0.04815074428053346, 0.0008930454810143837, 0.05318462260012811, 0.0003808294014608455;

  } else if (context.isRank(2)) { //Slave2

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.773301313815872e-05, -1.903469331632798e-06, 1.873284151703964e-05, 1.007255996407046e-06, 1.939089962749955e-05, -8.922690925673692e-07, 2.037857848281338e-05, -3.216215773080572e-07, 0, 0, 0, 0, 7.402516024117914e-06, 2.838268713863215e-06, 6.249761203424372e-06, -3.263999161516882e-06, 1.437336512054573e-05, 2.390776372597577e-06, 1.338687393172981e-05, -3.023290389839287e-06, 1.875511666670207e-05, -3.278415622262174e-07;
    fref << -0.01540107886229656, -0.0002615855206049247, 0.05487671246695029, 0.0002368365301820234, -0.0160990505051311, -0.0002741605822288765, 0.05415687969310452, 0.0002743355004921105, -0.007162227035058467, -0.007120294400987128, 0.02796795558583163, -0.02802370793272083, 0.05586613813214418, 0.0001029263016352541, -0.0144111353804976, -0.0001000512803033142, 0.05547743301530601, 0.0001393622825766435, -0.01479257484776803, -0.0001288175608001121, -0.01706063837893593, -0.0003583838471874983;
    drefNorm = 4.824205753272403e-05;
    frefNorm = 0.1225851627302816;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm), data.at(4)->values().norm());
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm), data.at(5)->values().norm());

    // update cplData
    displacements->values() << 1.716650969972045e-05, -1.856138836171773e-06, 1.818701485070425e-05, 9.439657883607802e-07, 1.874709534954619e-05, -8.85448704675396e-07, 1.975527973304359e-05, -3.501096287428596e-07, 0, 0, 0, 0, 7.228951427433641e-06, 2.745909101918556e-06, 6.052367643912141e-06, -3.179587143921995e-06, 1.398276918926419e-05, 2.29762824040882e-06, 1.297587398676e-05, -2.941551341709183e-06, 1.811863361465251e-05, -3.546317448342288e-07;
    doldValues << 1.773301313815872e-05, -1.903469331632798e-06, 1.873284151703964e-05, 1.007255996407046e-06, 1.939089962749955e-05, -8.922690925673692e-07, 2.037857848281338e-05, -3.216215773080572e-07, 0, 0, 0, 0, 7.402516024117914e-06, 2.838268713863215e-06, 6.249761203424372e-06, -3.263999161516882e-06, 1.437336512054573e-05, 2.390776372597577e-06, 1.338687393172981e-05, -3.023290389839287e-06, 1.875511666670207e-05, -3.278415622262174e-07;
    forces->values() << -0.09539527385890252, 0.0003208855941258066, -0.1853399184726223, 7.203155656644242e-05, -0.09532865072058605, 0.0009202649288056726, -0.1847925968312873, -0.0007589246108979722, -0.03998875591551594, -0.03982927597079221, -0.08489044889406808, 0.08470593806523596, -0.1739740974580442, 0.0007742373134568178, -0.08383286811708256, -0.0005911288917162662, -0.1811747642897668, 0.001020161732709184, -0.09112767929864005, -0.0008931566039992005, -0.08987332323372975, 0.002763113283891189;
    foldValues << -0.01540107886229656, -0.0002615855206049247, 0.05487671246695029, 0.0002368365301820234, -0.0160990505051311, -0.0002741605822288765, 0.05415687969310452, 0.0002743355004921105, -0.007162227035058467, -0.007120294400987128, 0.02796795558583163, -0.02802370793272083, 0.05586613813214418, 0.0001029263016352541, -0.0144111353804976, -0.0001000512803033142, 0.05547743301530601, 0.0001393622825766435, -0.01479257484776803, -0.0001288175608001121, -0.01706063837893593, -0.0003583838471874983;
  } else if (context.isRank(3)) { //Slave3
    // Dummy Slave to be able to reuse the 4 proc Master Slave Fixture
  }

  data.at(4)->oldValues.col(0) = doldValues;
  data.at(5)->oldValues.col(0) = foldValues;

  // QN- Update, 4. iteration
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master

  } else if (context.isRank(1)) { //Slave1

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.633368119919659e-05, -2.237360618564265e-07, 1.723960679736894e-05, -1.800451405634617e-06, 1.534489463502018e-05, -2.166297690952519e-06, 1.446268777631733e-05, -4.845152271528126e-07, 1.162685558347312e-05, -1.951619069146136e-06, 1.243725578988821e-05, -2.521440338918289e-06, 1.140405450734714e-05, -2.91813193235939e-06, 1.216867994540041e-05, -2.886191281742917e-06, 1.272486031358505e-05, -1.065380276784719e-06, 1.35816114080243e-05, -2.3781612091675e-06, 1.876166803225929e-05, -1.224865187618123e-06;
    fref << -0.01890549129941078, -0.0004389621840029381, 0.04954204963698766, 0.0004684469509178089, 0.04804708246519945, 0.0006388129115723739, -0.02036699235959898, -0.0007509963669595499, -0.02462352250264116, -0.00107121393884661, 0.04371708893453349, 0.001144538793507394, -0.01089395922122854, -0.01802447191842139, 0.02332094230979109, -0.01586002180275038, -0.02237733336181922, -0.0007488803918274927, 0.04617393036901896, 0.00081337820788483, 0.05067142953803216, 0.0003505856246817933;
    drefNorm = 4.815113092042934e-05;
    frefNorm = 0.1204042970796453;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm));
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm));

    // update cplData
    displacements->values() << 1.506845042291629e-05, -3.295713481574521e-07, 1.601708402767785e-05, -1.776032790440438e-06, 1.428998106709373e-05, -2.140925300298825e-06, 1.336604025915203e-05, -6.173668108734595e-07, 1.083614857997936e-05, -2.09020895349816e-06, 1.168515223592441e-05, -2.572621503366579e-06, 1.067901778881212e-05, -3.064421860106172e-06, 1.14804152344671e-05, -2.990689693425248e-06, 1.180274523671064e-05, -1.207361572630013e-06, 1.26994053163659e-05, -2.379420351559266e-06, 1.742665917249236e-05, -1.228726437307901e-06;
    doldValues << 1.633368119919659e-05, -2.237360618564265e-07, 1.723960679736894e-05, -1.800451405634617e-06, 1.534489463502018e-05, -2.166297690952519e-06, 1.446268777631733e-05, -4.845152271528126e-07, 1.162685558347312e-05, -1.951619069146136e-06, 1.243725578988821e-05, -2.521440338918289e-06, 1.140405450734714e-05, -2.91813193235939e-06, 1.216867994540041e-05, -2.886191281742917e-06, 1.272486031358505e-05, -1.065380276784719e-06, 1.35816114080243e-05, -2.3781612091675e-06, 1.876166803225929e-05, -1.224865187618123e-06;
    forces->values() << -0.07712271523301416, 0.004018032584755363, -0.1572746631951394, -0.003554957358717858, -0.1428154301136939, -0.005100471254411636, -0.06320679903752099, 0.00716455247252809, -0.02577587645775314, 0.009459639059138256, -0.1041793821597049, -0.008283348947413409, -0.02295652849878388, 0.02989440569403293, -0.06182269085435656, 0.01091909076145766, -0.04476922496853244, 0.00687112729333933, -0.1253800696382085, -0.006099818450316078, -0.1674291774005812, -0.002273881806091619;
    foldValues << -0.01890549129941078, -0.0004389621840029381, 0.04954204963698766, 0.0004684469509178089, 0.04804708246519945, 0.0006388129115723739, -0.02036699235959898, -0.0007509963669595499, -0.02462352250264116, -0.00107121393884661, 0.04371708893453349, 0.001144538793507394, -0.01089395922122854, -0.01802447191842139, 0.02332094230979109, -0.01586002180275038, -0.02237733336181922, -0.0007488803918274927, 0.04617393036901896, 0.00081337820788483, 0.05067142953803216, 0.0003505856246817933;

  } else if (context.isRank(2)) { //Slave2

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.686458723791834e-05, -1.810446627874213e-06, 1.781592254220799e-05, 9.575750096228968e-07, 1.844107125032693e-05, -8.488149913661512e-07, 1.938083795505493e-05, -3.062726369792662e-07, 0, 0, 0, 0, 7.040556208586022e-06, 2.699202297106055e-06, 5.943687851414432e-06, -3.104374082089934e-06, 1.367008353382051e-05, 2.273466384432095e-06, 1.27314174291305e-05, -2.875442657958392e-06, 1.783629704415144e-05, -3.12140240279707e-07;
    fref << -0.01626734826572977, -0.0002552779342950638, 0.05227538136126555, 0.0002350515292044683, -0.01695703999572442, -0.0002612258237118909, 0.05156927126324851, 0.0002631458821373055, -0.007517710152292339, -0.007474504607240594, 0.02674580053052369, -0.02680294718917284, 0.05337717466458169, 0.0001101957844225173, -0.01516291356738642, -0.0001053694865272746, 0.05291470170150131, 0.0001489003323108161, -0.01561921932892748, -0.0001370950096012615, -0.01784913805041409, -0.000324580522802835;
    drefNorm = 4.587992970636867e-05;
    frefNorm = 0.1180354804938519;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm));
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm));

    // update cplData
    displacements->values() << 1.563743909676446e-05, -1.707572586404205e-06, 1.663287551913161e-05, 8.210579991784308e-07, 1.704678071734513e-05, -8.336427145015805e-07, 1.803030552728031e-05, -3.673472962716038e-07, 0, 0, 0, 0, 6.663771864888832e-06, 2.499283366425937e-06, 5.516134032932667e-06, -2.921164340377279e-06, 1.282308279788757e-05, 2.07209067754735e-06, 1.184094543159743e-05, -2.698009821996337e-06, 1.645805878055576e-05, -3.696322259193852e-07;
    doldValues << 1.686458723791834e-05, -1.810446627874213e-06, 1.781592254220799e-05, 9.575750096228968e-07, 1.844107125032693e-05, -8.488149913661512e-07, 1.938083795505493e-05, -3.062726369792662e-07, 0, 0, 0, 0, 7.040556208586022e-06, 2.699202297106055e-06, 5.943687851414432e-06, -3.104374082089934e-06, 1.367008353382051e-05, 2.273466384432095e-06, 1.27314174291305e-05, -2.875442657958392e-06, 1.783629704415144e-05, -3.12140240279707e-07;
    forces->values() << -0.09143825207871489, 0.0002922798859936043, -0.1734744585823354, 8.018501629471556e-05, -0.09140909287296797, 0.0008614262538692869, -0.1729893406976169, -0.0007078492982043917, -0.03836482525057584, -0.03821118279703851, -0.07931609050916776, 0.07913795276507131, -0.1626218046186428, 0.0007411076261039719, -0.08039872451576649, -0.0005667402343291361, -0.1694857588798654, 0.0009766586358261331, -0.0873517838382746, -0.0008552683303008771, -0.08627064033821233, 0.002609015553424872;
    foldValues << -0.01626734826572977, -0.0002552779342950638, 0.05227538136126555, 0.0002350515292044683, -0.01695703999572442, -0.0002612258237118909, 0.05156927126324851, 0.0002631458821373055, -0.007517710152292339, -0.007474504607240594, 0.02674580053052369, -0.02680294718917284, 0.05337717466458169, 0.0001101957844225173, -0.01516291356738642, -0.0001053694865272746, 0.05291470170150131, 0.0001489003323108161, -0.01561921932892748, -0.0001370950096012615, -0.01784913805041409, -0.000324580522802835;
  } else if (context.isRank(3)) { //Slave3
    // Dummy Slave to be able to reuse the 4 proc Master Slave Fixture
  }

  data.at(4)->oldValues.col(0) = doldValues;
  data.at(5)->oldValues.col(0) = foldValues;

  // QN- Update, 5. iteration
  pp.performAcceleration(data);

  if (context.isMaster()) { //Master

  } else if (context.isRank(1)) { //Slave1

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.275776729912441e-06, -7.411719120649928e-07, 2.009844043916140e-06, -8.166362562036991e-07,
        1.98429549697312266297e-06, -1.006128466370864e-06, 1.26860076890459333335e-06, -1.038973060940335e-06,
        1.553926578960109e-06, -1.85501566670221e-06, 2.21293251735091021427e-06, -1.645357854619679e-06,
        1.820890655907848e-06, -2.41565468511731e-06, 2.44472211085250353034e-06, -2.153167028992053e-06,
        1.367377699262639e-06, -1.402360838056698e-06, 2.05877363493188315457e-06, -1.275601072931189e-06,
        2.16056223288899190422e-06, -6.758668464219906e-07;

    fref << -0.02494062387644205, 2.83457783442084623737e-05, 2.78926996928726099456e-02, 4.927873553950413e-05,
        0.02806450388830189, 4.133693821366594e-05, -0.02479702682061509, 7.828236584438153e-05,
        -0.02470268458079963, 3.17420733502629e-05, 0.0282235973672153, 0.0001624809599050109,
        -0.01213681562464092, -0.01299265397817983, 0.01440494173854578, -0.01303642911469127,
        -0.02467306876921166, 4.955128897708022e-05, 0.02820821342518329, 9.359794339403125e-05,
        0.02784148826297314, 7.700134509804057e-05;

    drefNorm = 7.910283453653413e-06;
    frefNorm = 0.08415800797163485;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i), 1e-8));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm));
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm));

  } else if (context.isRank(2)) { //Slave2

    dref = Eigen::VectorXd::Zero(22);
    fref = Eigen::VectorXd::Zero(22);

    dref << 1.782695896956317e-06, -3.609665456581163e-07, 2.549165865006206e-06, -2.9152246796606e-07,
        1.641256623136569e-06, -3.63485032734166e-07, 2.399541139942214e-06, -5.221960935703056e-07,
        0, 0, 0, 0,
        1.540562829986484e-06, 2.646900958343007e-07, 6.572995607065404e-07, -5.789378440843989e-07,
        2.312916074266845e-06, 2.909615579313369e-08, 1.55690286980020939237e-06, -4.907962897313283e-07,
        1.41306136610511245499e-06, -5.042593080795539e-07;

    fref << -0.02407925914468757, -0.0001962517462601011, 2.86388026412357221684e-02, 0.0002189694547629064,
        -0.02469057752137756, -0.0001420431374215673, 0.02806138636325711, 0.0001618659175134536,
        -0.01072215615486178, -0.01066770812537758, 0.01563815124261267, -0.01570783185165062,
        0.03075510899502418, 0.0001765574924817709, -0.02194129299700746, -0.0001523229943385313,
        0.02962481369149697, 2.35608194343154197722e-04, -0.02307508702308544, -0.0002109167490697018,
        -0.02495036462961023, -1.65224214831653608282e-05;

    drefNorm = 5.676684399367158e-06;
    frefNorm = 0.08353026170200345;

    // validate values
    for (int i = 0; i < data.at(4)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(4)->values()(i), dref(i)));

    for (int i = 0; i < data.at(5)->values().size(); i++)
      BOOST_TEST(testing::equals(data.at(5)->values()(i), fref(i)));

    // validate norm
    BOOST_TEST(testing::equals(data.at(4)->values().norm(), drefNorm));
    BOOST_TEST(testing::equals(data.at(5)->values().norm(), frefNorm));
  } else if (context.isRank(3)) { //Slave3
    // Dummy Slave to be able to reuse the 4 proc Master Slave Fixture
  }
}

/// Test that runs on 4 processors.
BOOST_AUTO_TEST_CASE(testColumnsLogging)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  double           initialRelaxation        = 0.1;
  int              maxIterationsUsed        = 3;
  int              timestepsReused          = 1;
  int              filter                   = BaseQNAcceleration::QR1FILTER;
  double           singularityLimit         = 0.1;
  bool             enforceInitialRelaxation = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  std::vector<double> factors;
  factors.resize(1.0, 1.0);
  PtrPreconditioner prec(new ConstantPreconditioner(factors));
  std::vector<int>  vertexOffsets{2, 3, 3, 4};

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, false, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  IQNILSAcceleration acc(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                         timestepsReused, filter, singularityLimit, dataIDs, prec);

  mesh::PtrData   displacements(new mesh::Data("dvalues", -1, 1));
  Eigen::VectorXd dcol1;

  DataMap data;

  if (context.isMaster()) { // 2 vertices
    //init displacements
    Eigen::VectorXd insert(2);
    insert << 1.0, 1.0;
    utils::append(displacements->values(), insert);
    insert << 0.5, 0.5;
    utils::append(dcol1, insert);
  } else if (context.isRank(1)) { //1 vertex
    Eigen::VectorXd insert(1);
    insert << 1.0;
    utils::append(displacements->values(), insert);
    insert << 0.5;
    utils::append(dcol1, insert);
  } else if (context.isRank(2)) { //no vertices
  } else {                        //1 vertex
    Eigen::VectorXd insert(1);
    insert << 1.0;
    utils::append(displacements->values(), insert);
    insert << 0.5;
    utils::append(dcol1, insert);
  }

  PtrCouplingData dpcd(new CouplingData(displacements, dummyMesh, false));
  data.insert(std::pair<int, PtrCouplingData>(0, dpcd));
  acc.initialize(data);
  dpcd->oldValues.col(0) = dcol1;

  acc.performAcceleration(data);

  Eigen::VectorXd newdvalues1;
  if (context.isMaster()) {
    utils::append(newdvalues1, 1.1);
    utils::append(newdvalues1, 1.0);
  } else if (context.isRank(1)) {
    utils::append(newdvalues1, 1.0);
  } else if (context.isRank(2)) {
  } else if (context.isRank(3)) {
    utils::append(newdvalues1, 1.0);
  }
  data.begin()->second->values() = newdvalues1;

  acc.performAcceleration(data);

  Eigen::VectorXd newdvalues2;
  if (context.isMaster()) {
    utils::append(newdvalues2, 1.0);
    utils::append(newdvalues2, 2.0);
  } else if (context.isRank(1)) {
    utils::append(newdvalues2, 1.0);
  } else if (context.isRank(2)) {
  } else if (context.isRank(3)) {
    utils::append(newdvalues2, 1.0);
  }
  data.begin()->second->values() = newdvalues2;

  acc.iterationsConverged(data);

  BOOST_TEST(acc.getLSSystemCols() == 2);
  BOOST_TEST(acc.getDeletedColumns() == 0);
  BOOST_TEST(acc.getDroppedColumns() == 0);

  Eigen::VectorXd newdvalues3;
  if (context.isMaster()) {
    utils::append(newdvalues3, 1.1);
    utils::append(newdvalues3, 2.0);
  } else if (context.isRank(1)) {
    utils::append(newdvalues3, 1.0);
  } else if (context.isRank(2)) {
  } else if (context.isRank(3)) {
    utils::append(newdvalues3, 1.0);
  }
  data.begin()->second->values() = newdvalues3;

  acc.performAcceleration(data);

  Eigen::VectorXd newdvalues4;
  if (context.isMaster()) {
    utils::append(newdvalues4, 1.0);
    utils::append(newdvalues4, 1.0);
  } else if (context.isRank(1)) {
    utils::append(newdvalues4, 1.0);
  } else if (context.isRank(2)) {
  } else if (context.isRank(3)) {
    utils::append(newdvalues4, 5.0);
  }
  data.begin()->second->values() = newdvalues4;

  acc.iterationsConverged(data);

  BOOST_TEST(acc.getLSSystemCols() == 1);
  BOOST_TEST(acc.getDeletedColumns() == 1);
  BOOST_TEST(acc.getDroppedColumns() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

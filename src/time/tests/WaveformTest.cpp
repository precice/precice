#include <Eigen/Core>
#include "com/MPIDirectCommunication.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Waveform.hpp"

using namespace precice;
using namespace precice::time;
using namespace precice::cplscheme;

BOOST_AUTO_TEST_SUITE(TimeTests)

BOOST_AUTO_TEST_SUITE(WaveformTests)

BOOST_AUTO_TEST_CASE(testExtrapolateData)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;

  PtrMesh mesh(new Mesh("MyMesh", 3, false, testing::nextMeshID()));
  PtrData data   = mesh->createData("MyData", 1);
  int     dataID = data->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  BOOST_TEST(data->values().size() == 1);

  double                maxTime      = CouplingScheme::UNDEFINED_TIME;
  int                   maxTimesteps = 1;
  double                dt           = 1.0;
  std::string           first        = "First";
  std::string           second       = "Second";
  std::string           accessor     = second;
  com::PtrCommunication com(new com::MPIDirectCommunication());
  m2n::PtrM2N           globalCom(new m2n::M2N(com, m2n::DistributedComFactory::SharedPointer()));
  int                   maxIterations = 1;

  // Test first order extrapolation
  SerialCouplingScheme scheme(maxTime, maxTimesteps, dt, 16, first, second,
                              accessor, globalCom, constants::FIXED_TIME_WINDOW_SIZE,
                              BaseCouplingScheme::Implicit, maxIterations);

  scheme.addDataToSend(data, mesh, true);
  scheme.setExtrapolationOrder(1);
  scheme.setupDataMatrices();
  Waveform *    waveform = scheme.getWaveform(dataID);
  CouplingData *cplData  = scheme.getSendData(dataID);
  BOOST_CHECK(waveform); // no nullptr
  BOOST_CHECK(cplData);  // no nullptr
  BOOST_TEST(cplData->values().size() == 1);
  BOOST_TEST(cplData->lastIteration().size() == 1);
  BOOST_TEST(waveform->lastTimeWindows().cols() == 2);
  BOOST_TEST(waveform->lastTimeWindows().rows() == 1);

  BOOST_TEST(testing::equals(cplData->values()(0), 0.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 0.0));

  cplData->values()(0) = 1.0;
  scheme.setTimeWindows(scheme.getTimeWindows() + 1);
  scheme.storeWindowData();
  scheme.extrapolateData();
  scheme.storeLastIteration();
  BOOST_TEST(testing::equals(cplData->values()(0), 2.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 2.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 0.0));

  cplData->values()(0) = 4.0;
  scheme.setTimeWindows(scheme.getTimeWindows() + 1);
  scheme.storeWindowData();
  scheme.extrapolateData();
  scheme.storeLastIteration();
  BOOST_TEST(testing::equals(cplData->values()(0), 7.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 7.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 1.0));

  // Test second order extrapolation
  cplData->values() = Eigen::VectorXd::Zero(cplData->values().size());
  cplData->storeIteration();
  SerialCouplingScheme scheme2(maxTime, maxTimesteps, dt, 16, first, second, accessor, globalCom, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, maxIterations);

  scheme2.addDataToSend(data, mesh, false);
  scheme2.setExtrapolationOrder(2);
  scheme2.setupDataMatrices();
  cplData  = scheme2.getSendData(dataID);
  waveform = scheme2.getWaveform(dataID);
  BOOST_CHECK(cplData); // no nullptr
  BOOST_TEST(cplData->values().size() == 1);
  BOOST_TEST(waveform->lastTimeWindows().cols() == 3);
  BOOST_TEST(waveform->lastTimeWindows().rows() == 1);
  BOOST_TEST(cplData->lastIteration().size() == 1);
  BOOST_TEST(testing::equals(cplData->values()(0), 0.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 2), 0.0));

  cplData->values()(0) = 1.0;
  scheme2.setTimeWindows(scheme2.getTimeWindows() + 1);
  scheme2.storeWindowData();
  scheme2.extrapolateData();
  scheme2.storeLastIteration();
  BOOST_TEST(testing::equals(cplData->values()(0), 2.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 2.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 2), 0.0));

  cplData->values()(0) = 4.0;
  scheme2.setTimeWindows(scheme2.getTimeWindows() + 1);
  scheme2.storeWindowData();
  scheme2.extrapolateData();
  scheme2.storeLastIteration();
  BOOST_TEST(testing::equals(cplData->values()(0), 8.0));
  BOOST_TEST(testing::equals(cplData->lastIteration()(0), 8.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform->lastTimeWindows()(0, 2), 0.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
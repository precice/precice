#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "xml/XMLTag.hpp"
#include <Eigen/Core>

#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

using namespace precice;
using namespace precice::cplscheme;

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

void runCoupling(
      CouplingScheme&                cplScheme,
      const std::string&             nameParticipant,
      const mesh::MeshConfiguration& meshConfig,
      const std::vector<int>&        validIterations )
{
  BOOST_TEST(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  BOOST_TEST(mesh->data().size() == 2);
  BOOST_TEST(mesh->vertices().size() > 0);
  mesh::Vertex& vertex = mesh->vertices()[0];
  int index = vertex.getID();
  auto& dataValues0 = mesh->data()[0]->values();
  auto& dataValues1 = mesh->data()[1]->values();
  double initialStepsizeData0 = 5.0;
  double stepsizeData0 = 5.0;
  Eigen::VectorXd initialStepsizeData1 = Eigen::VectorXd::Constant(3, 5.0);
  Eigen::VectorXd stepsizeData1 = Eigen::VectorXd::Constant(3, 5.0);
  double computedTime = 0.0;
  int computedTimesteps = 0;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  BOOST_TEST (((nameParticipant == nameParticipant0) || (nameParticipant == nameParticipant1)));
  int iterationCount = 0;
  std::vector<int>::const_iterator iterValidIterations = validIterations.begin();

  if( nameParticipant == nameParticipant0 ){
    cplScheme.initialize ( 0.0, 1);
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(not cplScheme.hasDataBeenExchanged());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction(constants::actionWriteIterationCheckpoint());
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));

    while ( cplScheme.isCouplingOngoing() ) {
      dataValues0[index] += stepsizeData0;
      // The max timestep length is required to be obeyed.
      double maxLengthTimestep = cplScheme.getNextTimestepMaxLength();
      cplScheme.addComputedTime ( maxLengthTimestep);
      cplScheme.advance();
      iterationCount++;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ) {
        // Advance participant time and timestep
        computedTime += maxLengthTimestep;
        computedTimesteps ++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime() ));
        BOOST_TEST(testing::equals(computedTimesteps, cplScheme.getTimesteps()-1 ));
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(testing::equals(iterationCount, *iterValidIterations ));
        if ( cplScheme.isCouplingOngoing() ) {
          BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction(constants::actionWriteIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
        }
        else {
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint() ));
        }
        iterationCount = 0;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          BOOST_TEST(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData0 = initialStepsizeData0;
      }
      else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        BOOST_TEST(iterationCount < *iterValidIterations);
        BOOST_TEST(cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
        cplScheme.performedAction ( constants::actionReadIterationCheckpoint());
        BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        stepsizeData0 -= 1.0;
      }
      // the first participant always receives new data
      //if(cplScheme.isCouplingOngoing())
      BOOST_TEST(cplScheme.hasDataBeenExchanged());
    }
    cplScheme.finalize (); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3 ));
    BOOST_TEST(testing::equals(computedTimesteps, 3 ));
  }
  else if (nameParticipant == nameParticipant1) {
    cplScheme.initialize ( 0.0, 1);
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(cplScheme.hasDataBeenExchanged());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( constants::actionWriteIterationCheckpoint());
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));

    while ( cplScheme.isCouplingOngoing() ) {
      Eigen::VectorXd currentData(3);
      currentData = dataValues1.segment(index * 3, 3);
      currentData += stepsizeData1;
      dataValues1.segment(index * 3, 3) = currentData;
      // The max timestep length is required to be obeyed.
      double maxLengthTimestep = cplScheme.getNextTimestepMaxLength();
      cplScheme.addComputedTime ( maxLengthTimestep);
      cplScheme.advance();
      iterationCount++;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ) {
        // Advance participant time and timestep
        computedTime += maxLengthTimestep;
        computedTimesteps ++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime() ));
        BOOST_TEST(testing::equals(computedTimesteps, cplScheme.getTimesteps()-1 ));
        // The iterations are enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(testing::equals(iterationCount, *iterValidIterations ));
        if ( cplScheme.isCouplingOngoing() ) {
          BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction ( constants::actionWriteIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
        }
        else {
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint() ));
        }
        iterationCount = 0;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          BOOST_TEST(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData1 = initialStepsizeData1;
      }
      else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        BOOST_TEST(iterationCount < *iterValidIterations);
        BOOST_TEST(cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
        // The load checkpoint action requires to fallback to the cplScheme of the
        // first implicit iteration of the current timestep/time.
        cplScheme.performedAction ( constants::actionReadIterationCheckpoint());
        BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        //stepsizeData1 -= 1.0;
        stepsizeData1 -= Eigen::Vector3d::Constant(1.0);
      }
      // only check if data is received
      if(cplScheme.isCouplingOngoing())
        BOOST_TEST(cplScheme.hasDataBeenExchanged());
    }
    cplScheme.finalize (); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3 ));
    BOOST_TEST(testing::equals(computedTimesteps, 3 ));
  }
}

void runCouplingWithSubcycling
(
    CouplingScheme&                cplScheme,
    const std::string&             nameParticipant,
    const mesh::MeshConfiguration& meshConfig,
    const std::vector<int>&        validIterations )
{
  BOOST_TEST(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  BOOST_TEST(mesh->data().size() == 2);
  BOOST_TEST(mesh->vertices().size() > 0);
  double initialStepsizeData0 = 5.0;
  double stepsizeData0 = 5.0;
  Eigen::Vector3d initialStepsizeData1 = Eigen::Vector3d::Constant ( 5.0);
  Eigen::Vector3d stepsizeData1 = Eigen::Vector3d::Constant( 5.0);
  double computedTime = 0.0;
  int computedTimesteps = 0;
  std::string nameParticipant0 ( "Participant0");
  std::string nameParticipant1 ( "Participant1");
  BOOST_TEST(((nameParticipant == nameParticipant0) || (nameParticipant == nameParticipant1)));
  int iterationCount = 0;
  std::vector<int>::const_iterator iterValidIterations =
      validIterations.begin();

  if ( nameParticipant == nameParticipant0 ) {
    iterationCount++; // different handling due to subcycling
    cplScheme.initialize ( 0.0, 1);
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(not cplScheme.hasDataBeenExchanged());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( constants::actionWriteIterationCheckpoint());
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));

    double maxTimestepLength = cplScheme.getNextTimestepMaxLength();
    double computedTimestepLength = maxTimestepLength / 2.0;
    int subcyclingStep = 0;

    // Main coupling loop
    while ( cplScheme.isCouplingOngoing() ){
      cplScheme.addComputedTime ( computedTimestepLength);
      cplScheme.advance();
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ){
        // Advance participant time and timestep
        computedTime += maxTimestepLength;
        computedTimesteps ++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(testing::equals(computedTimesteps, cplScheme.getTimesteps()-1));
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(testing::equals(iterationCount, *iterValidIterations));
        if ( cplScheme.isCouplingOngoing() ) {
          BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction(constants::actionWriteIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
        }
        else {
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint() ));
        }
        iterationCount = 1;
        iterValidIterations++;
        if (iterValidIterations == validIterations.end()) {
          BOOST_TEST(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData0 = initialStepsizeData0;
        BOOST_TEST(testing::equals(subcyclingStep, 1));
        subcyclingStep = 0;
      }
      else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        // If length of global timestep is reached
        if ( cplScheme.hasDataBeenExchanged() ) {
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction ( constants::actionReadIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData0 -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++; // Implicit coupling iterations
        }
        else { // If subcycling
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(subcyclingStep < 2);
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize (); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(testing::equals(computedTimesteps, 3));
  }

  else if ( nameParticipant == nameParticipant1 ) {
    iterationCount++; // different handling due to subcycling
    cplScheme.initialize ( 0.0, 1);
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(cplScheme.hasDataBeenExchanged());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( constants::actionWriteIterationCheckpoint());
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));

    double maxTimestepLength = cplScheme.getNextTimestepMaxLength ();
    double preferredTimestepLength = maxTimestepLength / 2.5;
    double computedTimestepLength = preferredTimestepLength;
    int subcyclingStep = 0;

    // Main coupling loop
    while ( cplScheme.isCouplingOngoing() ){
      cplScheme.addComputedTime ( computedTimestepLength);
      cplScheme.advance();
      computedTimestepLength =
          cplScheme.getNextTimestepMaxLength() < preferredTimestepLength
          ? cplScheme.getNextTimestepMaxLength()
              : preferredTimestepLength;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ){
        // Advance participant time and timestep
        computedTime += maxTimestepLength;
        computedTimesteps ++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(testing::equals(computedTimesteps, cplScheme.getTimesteps()-1));
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(testing::equals(iterationCount, *iterValidIterations));
        if ( cplScheme.isCouplingOngoing() ) {
          BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction(constants::actionWriteIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
        }
        else {
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint() ));
        }
        iterationCount = 1;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          BOOST_TEST(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData1 = initialStepsizeData1;
        BOOST_TEST(testing::equals(subcyclingStep, 2));
        subcyclingStep = 0;
      }
      else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        // If length of global timestep is reached
        if ( cplScheme.hasDataBeenExchanged() ) {
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          cplScheme.performedAction ( constants::actionReadIterationCheckpoint());
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData1.array() -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++; // Implicit coupling iterations
        }
        else { // If subcycling
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
          BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
          BOOST_TEST(subcyclingStep < 3);
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize (); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(testing::equals(computedTimesteps, 3));
  }
}

struct SerialImplicitCouplingSchemeFixture
{
  std::string _pathToTests;

  SerialImplicitCouplingSchemeFixture(){
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }
};

BOOST_FIXTURE_TEST_SUITE(SerialImplicitCouplingSchemeTests, SerialImplicitCouplingSchemeFixture)

BOOST_AUTO_TEST_CASE(testParseConfigurationWithRelaxation)
{
  using namespace mesh;

  std::string path(_pathToTests + "serial-implicit-cplscheme-relax-const-config.xml");

  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, xml::ConfigurationContext{}, path);
  BOOST_CHECK(cplSchemeConfig._accelerationConfig->getAcceleration().get()); // no nullptr
  meshConfig->setMeshSubIDs();
}

BOOST_AUTO_TEST_CASE(testExtrapolateData)
{
  using namespace mesh;

  PtrMesh mesh(new Mesh("MyMesh", 3, false));
  PtrData data = mesh->createData("MyData", 1);
  int dataID = data->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  BOOST_TEST(data->values().size() == 1);

  double maxTime = CouplingScheme::UNDEFINED_TIME;
  int maxTimesteps = 1;
  double dt = 1.0;
  std::string first = "First";
  std::string second = "Second";
  std::string accessor = second;
  com::PtrCommunication com(new com::MPIDirectCommunication());
  m2n::PtrM2N globalCom(new m2n::M2N(com, m2n::DistributedComFactory::SharedPointer()));
  int maxIterations = 1;

  // Test first order extrapolation
  SerialCouplingScheme scheme(maxTime, maxTimesteps, dt, 16, first, second,
      accessor, globalCom, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, maxIterations);

  scheme.addDataToSend(data, mesh, true);
  scheme.setExtrapolationOrder(1);
  scheme.setupDataMatrices(scheme.getSendData());
  CouplingData* cplData = scheme.getSendData(dataID);
  BOOST_CHECK(cplData); // no nullptr
  BOOST_TEST(cplData->values->size() == 1);
  BOOST_TEST(cplData->oldValues.cols() == 2);
  BOOST_TEST(cplData->oldValues.rows() == 1);
  BOOST_TEST(testing::equals((*cplData->values)[0], 0.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 0.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 0.0));

  (*cplData->values)[0] = 1.0;
  scheme.setTimesteps(scheme.getTimesteps() + 1);
  scheme.extrapolateData(scheme.getSendData());
  BOOST_TEST(testing::equals((*cplData->values)[0], 2.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 2.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 1.0));

  (*cplData->values)[0] = 4.0;
  scheme.setTimesteps(scheme.getTimesteps() + 1);
  scheme.extrapolateData(scheme.getSendData());
  BOOST_TEST(testing::equals((*cplData->values)[0], 7.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 7.0));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 4.0));

  // Test second order extrapolation
  *cplData->values = Eigen::VectorXd::Zero(cplData->values->size());
  cplData->oldValues = Eigen::MatrixXd::Zero(cplData->oldValues.rows(), cplData->oldValues.cols());
  //assign(*cplData->values) = 0.0;
  //assign(cplData->oldValues) = 0.0;
  SerialCouplingScheme scheme2(maxTime, maxTimesteps, dt, 16, first, second, accessor, globalCom, constants::FIXED_DT, BaseCouplingScheme::Implicit, maxIterations);

  scheme2.addDataToSend(data, mesh, false);
  scheme2.setExtrapolationOrder(2);
  scheme2.setupDataMatrices(scheme2.getSendData());
  cplData = scheme2.getSendData(dataID);
  BOOST_CHECK(cplData); // no nullptr
  BOOST_TEST(cplData->values->size() == 1);
  BOOST_TEST(cplData->oldValues.cols() == 3);
  BOOST_TEST(cplData->oldValues.rows() == 1);
  BOOST_TEST(testing::equals((*cplData->values)[0], 0.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 0.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 0.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,2), 0.0 ));

  (*cplData->values)[0] = 1.0;
  scheme2.setTimesteps(scheme2.getTimesteps() + 1);
  scheme2.extrapolateData(scheme2.getSendData());
  BOOST_TEST(testing::equals((*cplData->values)[0], 2.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 2.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 1.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,2), 0.0 ));

  (*cplData->values)[0] = 4.0;
  scheme2.setTimesteps(scheme2.getTimesteps() + 1);
  scheme2.extrapolateData(scheme2.getSendData());
  BOOST_TEST(testing::equals((*cplData->values)[0], 8.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,0), 8.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,1), 4.0 ));
  BOOST_TEST(testing::equals(cplData->oldValues(0,2), 1.0 ));
}

/// Test that runs on 2 processors.
BOOST_FIXTURE_TEST_CASE(testAbsConvergenceMeasureSynchronized, testing::M2NFixture,
                      * testing::MinRanks(2)
                      * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData("data0", 1);
  dataConfig->addData("data1", 3);

  MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new Mesh("Mesh", 3, false));
  mesh->createData("data0", 1);
  mesh->createData("data1", 3);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues ();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  double maxTime = 1.0;
  int maxTimesteps = 3;
  double timestepLength = 0.1;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  std::string nameLocalParticipant("");
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  if (utils::Parallel::getProcessRank() == 0) {
    nameLocalParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
  }
  else if (utils::Parallel::getProcessRank() == 1 ) {
    nameLocalParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
  }

  // Create the coupling scheme object
  cplscheme::SerialCouplingScheme cplScheme (
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0,
      nameParticipant1, nameLocalParticipant, m2n, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, 100
  );
  cplScheme.addDataToSend (mesh->data()[sendDataIndex], mesh, false);
  cplScheme.addDataToReceive (mesh->data()[receiveDataIndex], mesh, false);

  double convergenceLimit1 = sqrt(3.0); // when diff_vector = (1.0, 1.0, 1.0)
  cplscheme::impl::PtrConvergenceMeasure absoluteConvMeasure1 (
      new cplscheme::impl::AbsoluteConvergenceMeasure(convergenceLimit1));
  cplScheme.addConvergenceMeasure(mesh->data()[1], false, false, absoluteConvMeasure1);

  // Expected iterations per implicit timesptep
  std::vector<int> validIterations = {5, 5, 5};
  runCoupling (cplScheme, nameLocalParticipant, meshConfig, validIterations);
}

BOOST_AUTO_TEST_CASE(testConfiguredAbsConvergenceMeasureSynchronized,
                   * testing::MinRanks(2)
                   * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;

  std::string configurationPath (
      _pathToTests + "serial-implicit-cplscheme-absolute-config.xml");

  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, xml::ConfigurationContext{}, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("Participant0", "Participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  std::vector<int> validIterations = {5, 5, 5};

  std::string nameLocalParticipant ("");
  if (utils::Parallel::getProcessRank() == 0 ) {
    nameLocalParticipant = "Participant0";
    utils::Parallel::splitCommunicator(nameLocalParticipant);
    m2n->requestMasterConnection ("Participant1", "Participant0");
  }
  else if (utils::Parallel::getProcessRank() == 1 ) {
    nameLocalParticipant = "Participant1";
    utils::Parallel::splitCommunicator(nameLocalParticipant);
    m2n->acceptMasterConnection ("Participant1", "Participant0");
  }

  runCoupling (*cplSchemeConfig.getCouplingScheme(nameLocalParticipant),
      nameLocalParticipant, *meshConfig, validIterations);
  utils::Parallel::clearGroups();
}

BOOST_FIXTURE_TEST_CASE(testMinIterConvergenceMeasureSynchronized, testing::M2NFixture,
                      * testing::MinRanks(2)
                      * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig (new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData ("data0", 1);
  dataConfig->addData ("data1", 3);

  mesh::MeshConfiguration meshConfig (root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh (new mesh::Mesh("Mesh", 3, false));
  mesh->createData ("data0", 1);
  mesh->createData ("data1", 3);
  mesh->createVertex (Eigen::Vector3d::Zero());
  mesh->allocateDataValues ();
  meshConfig.addMesh (mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  double maxTime = 1.0;
  int maxTimesteps = 3;
  double timestepLength = 0.1;
  std::string nameParticipant0 ("Participant0");
  std::string nameParticipant1 ("Participant1");
  std::string nameLocalParticipant ("");
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  if (utils::Parallel::getProcessRank() == 0 ) {
    nameLocalParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
  }
  else if (utils::Parallel::getProcessRank() == 1 ) {
    nameLocalParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
  }

  // Create the coupling scheme object
  cplscheme::SerialCouplingScheme cplScheme (
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
      nameLocalParticipant, m2n, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend (mesh->data()[sendDataIndex], mesh, false);
  cplScheme.addDataToReceive (mesh->data()[receiveDataIndex], mesh, false);

  // Add convergence measures
  int minIterations = 3;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1 (
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(mesh->data()[1], false, false, minIterationConvMeasure1);

  // Expected iterations per implicit timesptep
  std::vector<int> validIterations = {3, 3, 3};
  runCoupling (cplScheme, nameLocalParticipant, meshConfig, validIterations);
}

BOOST_FIXTURE_TEST_CASE(testMinIterConvergenceMeasureSynchronizedWithSubcycling, testing::M2NFixture,
                      * testing::MinRanks(2)
                      * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData ( "data0", 1);
  dataConfig->addData ( "data1", 3);

  mesh::MeshConfiguration meshConfig ( root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh ( new mesh::Mesh("Mesh", 3, false));
  mesh->createData ( "data0", 1);
  mesh->createData ( "data1", 3);
  mesh->createVertex ( Eigen::Vector3d::Zero());
  mesh->allocateDataValues ();
  meshConfig.addMesh ( mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  double maxTime = 1.0;
  int maxTimesteps = 3;
  double timestepLength = 0.1;
  std::string nameParticipant0 ( "Participant0");
  std::string nameParticipant1 ( "Participant1");
  std::string nameLocalParticipant ( "");
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  std::vector<int> validIterations;
  if ( utils::Parallel::getProcessRank() == 0 ) {
    nameLocalParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
    validIterations = {3, 3, 3};
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    nameLocalParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
    validIterations = {3, 3, 3};
  }

  // Create the coupling scheme object
  cplscheme::SerialCouplingScheme cplScheme (
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
      nameLocalParticipant, m2n, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh, false);
  cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh, false);

  // Add convergence measures
  int minIterations = 3;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1 (
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(mesh->data()[1], false, false, minIterationConvMeasure1);
  runCouplingWithSubcycling (
      cplScheme, nameLocalParticipant, meshConfig, validIterations);
}

BOOST_FIXTURE_TEST_CASE(testInitializeData, testing::M2NFixture,
                   * testing::MinRanks(2)
                   * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  xml::XMLTag root = xml::getRootTag();

  // Create a data configuration, to simplify configuration of data

  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
  const auto dataID0 = mesh->createData("Data0", 1)->getID();
  const auto dataID1 = mesh->createData("Data1", 3)->getID();;
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  double maxTime = 1.0;
  int maxTimesteps = 3;
  double timestepLength = 0.1;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  std::string nameLocalParticipant("");
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  bool initData = false;
  if (utils::Parallel::getProcessRank() == 0){
    nameLocalParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
  }
  else if (utils::Parallel::getProcessRank() == 1){
    nameLocalParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
    initData = true;
  }

  // Create the coupling scheme object
  cplscheme::SerialCouplingScheme cplScheme(
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
      nameLocalParticipant, m2n, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend(mesh->data()[sendDataIndex], mesh, initData);
  cplScheme.addDataToReceive(mesh->data()[receiveDataIndex], mesh, not initData);

  // Add convergence measures
  int minIterations = 3;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1 (
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(mesh->data()[1], false, false, minIterationConvMeasure1);

  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  cplScheme.initialize(0.0, 1);

  if (nameLocalParticipant == nameParticipant0){
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenExchanged());
    auto& values = mesh->data(dataID1)->values();
    BOOST_TEST(testing::equals(values, Eigen::Vector3d(1.0, 2.0, 3.0)));
    mesh->data(dataID0)->values() = Eigen::VectorXd::Constant(mesh->data(dataID0)->values().size(), 4.0);
    while (cplScheme.isCouplingOngoing()){
      if (cplScheme.isActionRequired(writeIterationCheckpoint)){
        cplScheme.performedAction(writeIterationCheckpoint);
      }
      if (cplScheme.isActionRequired(readIterationCheckpoint)){
        cplScheme.performedAction(readIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
    }
  }
  else {
    BOOST_TEST(nameLocalParticipant == nameParticipant1);
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    cplScheme.performedAction(constants::actionWriteInitialData());
    auto& values = mesh->data(dataID0)->values();
    BOOST_TEST(testing::equals(values(0), 0.0));
    Eigen::VectorXd v(3); v << 1.0, 2.0, 3.0;
    mesh->data(dataID1)->values() = v;
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenExchanged());
    BOOST_TEST(testing::equals(values(0), 4.0));
    while (cplScheme.isCouplingOngoing()){
      if (cplScheme.isActionRequired(writeIterationCheckpoint)){
        cplScheme.performedAction(writeIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
      if (cplScheme.isActionRequired(readIterationCheckpoint)){
        cplScheme.performedAction(readIterationCheckpoint);
      }
    }
  }
  cplScheme.finalize();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

# endif // not PRECICE_NO_MPI

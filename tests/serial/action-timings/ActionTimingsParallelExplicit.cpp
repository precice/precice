#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <action/RecorderAction.hpp>
#include <precice/precice.hpp>
#include <vector>

/**
 * @brief Test to make sure that actions are called in the right order for explicit coupling via RecorderAction
 */
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(ActionTimings)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ActionTimingsParallelExplicit)
{
  PRECICE_TEST();

  using namespace precice;
  Participant interface(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  double      writeValue;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Forces";
    readDataName  = "Velocities";
    writeValue    = 1;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Velocities";
    readDataName  = "Forces";
    writeValue    = 2;
  }
  int                 dimensions = interface.getMeshDimensions(meshName);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = interface.setMeshVertex(meshName, vertex);

  double dt = -1;
  BOOST_TEST(action::RecorderAction::records.empty());
  action::RecorderAction::reset();
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);

  if (interface.requiresInitialData()) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    interface.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
  }

  interface.initialize();
  dt = interface.getMaxTimeStepSize();
  BOOST_TEST(dt == 1.0);
  BOOST_TEST(action::RecorderAction::records.size() == 2);
  BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_POST);
  BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::READ_MAPPING_POST);
  action::RecorderAction::reset();

  while (interface.isCouplingOngoing()) {
    interface.readData(meshName, readDataName, {&vertexID, 1}, dt, readData);
    interface.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
    interface.advance(dt);
    BOOST_TEST(interface.isTimeWindowComplete());

    BOOST_TEST(action::RecorderAction::records.size() == 2);
    BOOST_TEST(action::RecorderAction::records.at(0).timing == action::Action::WRITE_MAPPING_POST);
    BOOST_TEST(action::RecorderAction::records.at(1).timing == action::Action::READ_MAPPING_POST);

    action::RecorderAction::reset();
  }
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ActionTimings

#endif // PRECICE_NO_MPI

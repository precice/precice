#ifndef PRECICE_NO_MPI

#include "precice/impl/ParticipantImpl.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Whitebox)
/**
 * @brief Runs a coupled sim. with data scaling applied.
 *
 * SolverOne writes vector data on a cube geometry. The data values are defined
 * and stay constant over the coupling cycles. SolverTwo has a scaling of the
 * values activated and reads the scaled values.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(TestExplicitWithDataScaling)
{
  PRECICE_TEST();

  Participant cplInterface(context.name, context.config(), 0, 1);

  std::vector<double> positions = {0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.0};
  std::vector<int>    ids       = {0, 0, 0, 0};

  if (context.isNamed("SolverOne")) {
    auto meshName = "Test-Square-One";
    BOOST_REQUIRE(cplInterface.getMeshDimensions(meshName));
    cplInterface.setMeshVertices(meshName, positions, ids);
    for (int i = 0; i < 4; i++)
      cplInterface.setMeshEdge(meshName, ids.at(i), ids.at((i + 1) % 4));

    cplInterface.initialize();
    double dt = cplInterface.getMaxTimeStepSize();

    auto velocitiesID = "Velocities";
    while (cplInterface.isCouplingOngoing()) {
      std::vector<double> data;
      for (size_t i = 0; i < testing::WhiteboxAccessor::impl(cplInterface).mesh("Test-Square-One").nVertices(); ++i) {
        data.push_back(i);
        data.push_back(i);
      }
      cplInterface.writeData(meshName, velocitiesID, ids, data);
      cplInterface.advance(dt);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshName = "Test-Square-Two";
    BOOST_REQUIRE(cplInterface.getMeshDimensions(meshName));
    cplInterface.setMeshVertices(meshName, positions, ids);
    for (int i = 0; i < 4; i++)
      cplInterface.setMeshEdge(meshName, ids.at(i), ids.at((i + 1) % 4));

    cplInterface.initialize();
    double dt = cplInterface.getMaxTimeStepSize();

    auto velocitiesID = "Velocities";
    while (cplInterface.isCouplingOngoing()) {
      const auto size = testing::WhiteboxAccessor::impl(cplInterface).mesh("Test-Square-Two").nVertices();

      std::vector<double> expected;
      for (size_t i = 0; i < size; ++i) {
        expected.push_back(i * 10.0);
        expected.push_back(i * 10.0);
      }
      std::vector<double> data(expected.size());
      cplInterface.readData(meshName, velocitiesID, ids, dt, data);
      BOOST_TEST(data == expected, boost::test_tools::per_element());

      cplInterface.advance(dt);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Whitebox

#endif // PRECICE_NO_MPI

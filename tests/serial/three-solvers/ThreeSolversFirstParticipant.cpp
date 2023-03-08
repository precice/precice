#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(ThreeSolvers)
BOOST_AUTO_TEST_CASE(ThreeSolversFirstParticipant)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));
  std::string config = context.config();

  // SolverOne prescribes these, thus SolverTwo and SolverThree expect these (we use "first-participant" as dt method)
  std::vector<double> timestepSizes{1.0, 2.0, 3.0};

  precice::SolverInterface precice(context.name, config, 0, 1);

  if (context.isNamed("SolverOne")) {

    auto meshName = "MeshOne";
    precice.setMeshVertex(meshName, Eigen::Vector2d(0, 0).data());

    precice.initialize();

    for (auto dt : timestepSizes) {
      BOOST_TEST(precice.isCouplingOngoing());
      precice.advance(dt);
    }
    BOOST_TEST(not precice.isCouplingOngoing());
    precice.finalize();

  } else if (context.isNamed("SolverTwo")) {

    auto meshAID = "MeshTwoA";
    precice.setMeshVertex(meshAID, Eigen::Vector2d(0, 0).data());
    auto meshBID = "MeshTwoB";
    precice.setMeshVertex(meshBID, Eigen::Vector2d(0, 0).data());

    double dt = precice.initialize();

    for (auto expected_dt : timestepSizes) {
      BOOST_TEST(precice.isCouplingOngoing());
      BOOST_TEST(dt == expected_dt);
      dt = precice.advance(dt);
    }

    BOOST_TEST(not precice.isCouplingOngoing());
    precice.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverThree"));

    auto meshName = "MeshThree";
    precice.setMeshVertex(meshName, Eigen::Vector2d(0, 0).data());

    double dt = precice.initialize();

    for (auto expected_dt : timestepSizes) {
      BOOST_TEST(precice.isCouplingOngoing());
      BOOST_TEST(dt == expected_dt);
      dt = precice.advance(dt);
    }
    BOOST_TEST(not precice.isCouplingOngoing());
    precice.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ThreeSolvers

#endif // PRECICE_NO_MPI

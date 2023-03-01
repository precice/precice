#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/SolverInterface.hpp"

/**
 * @brief helper function for a simple test with data initialization
 */
void testDataInitialization(precice::testing::TestContext context, std::string config)
{
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, config, 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshID = "MeshOne";
    Vector3d pos    = Vector3d::Zero();
    cplInterface.setMeshVertex(meshID, pos.data());
    auto   dataID     = "Data"; //  meshOneID
    double valueDataB = 0.0;
    cplInterface.initialize();
    cplInterface.readScalarData(meshID, dataID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshID = "MeshTwo";
    Vector3d pos    = Vector3d::Zero();
    cplInterface.setMeshVertex(meshID, pos.data());

    //tell preCICE that data has been written
    BOOST_REQUIRE(cplInterface.requiresInitialData());

    auto dataID = "Data"; //  meshTwoID
    cplInterface.writeScalarData(meshID, dataID, 0, 2.0);
    cplInterface.initialize();
    cplInterface.finalize();
  }
}

#endif

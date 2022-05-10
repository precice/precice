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
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshOneID, pos.data());
    cplInterface.initialize();
    int    dataID     = cplInterface.getDataID("Data", meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshTwoID, pos.data());
    cplInterface.initialize();
    int dataID = cplInterface.getDataID("Data", meshTwoID);
    cplInterface.writeScalarData(dataID, 0, 2.0);
    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    cplInterface.finalize();
  }
}

#endif

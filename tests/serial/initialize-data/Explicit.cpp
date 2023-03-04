#include <boost/test/tools/old/interface.hpp>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(InitializeData)

/**
 * @brief The second solver initializes the data of the first.
 *
 * A mapping is employed for the second solver, i.e., at the end of
 * initialize(), the mapping needs to be invoked.
 */
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    cplInterface.setMeshVertex(meshName, Vector3d(1.0, 2.0, 3.0).data());
    auto   dataAID    = "DataOne"; //  meshOneID
    auto   dataBID    = "DataTwo"; //  meshOneID
    double valueDataB = 0.0;
    double maxDt      = cplInterface.initialize();
    cplInterface.readScalarData(meshName, dataBID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(meshName, dataAID, 0, valueDataA.data());
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readScalarData(meshName, dataBID, 0, valueDataB);
      BOOST_TEST(2.5 == valueDataB);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshName = "MeshTwo";
    Vector3d pos      = Vector3d::Zero();
    cplInterface.setMeshVertex(meshName, pos.data());

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    auto dataAID = "DataOne"; //  meshTwoID
    auto dataBID = "DataTwo"; //  meshTwoID
    cplInterface.writeScalarData(meshName, dataBID, 0, 2.0);
    //tell preCICE that data has been written and call initializeData
    double   maxDt = cplInterface.initialize();
    Vector3d valueDataA;
    cplInterface.readVectorData(meshName, dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName, dataBID, 0, 2.5);
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(meshName, dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // InitializeData

#endif // PRECICE_NO_MPI

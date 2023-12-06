#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)
BOOST_AUTO_TEST_CASE(RadialGeoMultiscale)
{
  PRECICE_TEST("Fluid1D"_on(1_rank), "Solid3D"_on(1_rank));
  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);
  if (context.isNamed("Fluid1D")) {
    Vector3d            coordOneA{0.0, 0.0, 0.0};
    Vector3d            coordOneB{0.0, 0.0, 1.0};
    Vector3d            coordOneC{0.0, 0.0, 2.0};
    std::vector<double> values;
    const unsigned int  nCoords = 3;
    for (unsigned int i = 0; i < nCoords; ++i) {
      values.emplace_back(std::pow(i + 1, 2));
    }
    auto             meshOneName = "Mesh1D";
    std::vector<int> ids;
    ids.emplace_back(cplInterface.setMeshVertex(meshOneName, coordOneA));
    ids.emplace_back(cplInterface.setMeshVertex(meshOneName, coordOneB));
    ids.emplace_back(cplInterface.setMeshVertex(meshOneName, coordOneC));

    auto dataAName = "HeatFluxLike";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    BOOST_TEST(cplInterface.isCouplingOngoing());

    cplInterface.writeData(meshOneName, dataAName, ids, values);
    cplInterface.advance(maxDt);

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Solid3D"));
    auto meshTwoName = "Mesh3D";

    Vector3d coordTwoA{0.5, 2.0, 0.0};
    Vector3d coordTwoB{0.5, 3.0, 1.0};
    Vector3d coordTwoC{0.5, 4.0, 2.0};

    // Setup receiving mesh.
    int idA = cplInterface.setMeshVertex(meshTwoName, coordTwoA);
    int idB = cplInterface.setMeshVertex(meshTwoName, coordTwoB);
    int idC = cplInterface.setMeshVertex(meshTwoName, coordTwoC);

    auto dataAName = "HeatFluxLike";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    double values[3];
    int    ids[] = {idA, idB, idC};

    cplInterface.readData(meshTwoName, dataAName, ids, maxDt, values);
    cplInterface.advance(maxDt);

    BOOST_TEST(values[0] == 1);
    BOOST_TEST(values[1] == 4);
    BOOST_TEST(values[2] == 9);

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // GeomultiscaleMapping

#endif // PRECICE_NO_MPI

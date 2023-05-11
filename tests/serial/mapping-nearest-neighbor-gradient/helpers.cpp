#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/Participant.hpp"
#include "testing/Testing.hpp"

using namespace precice;

void testVectorGradientFunctions(const TestContext &context)
{
  using Eigen::Vector3d;

  Participant interface(context.name, context.config(), 0, 1);
  if (context.isNamed("A")) {

    auto meshName = "MeshA";
    auto dataName = "DataA";

    Vector3d posOne = Vector3d::Constant(0.0);
    Vector3d posTwo = Vector3d::Constant(1.0);
    interface.setMeshVertex(meshName, posOne);
    interface.setMeshVertex(meshName, posTwo);

    // Initialize, thus sending the mesh.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double values[6]  = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int    indices[2] = {0, 1};
    interface.writeData(meshName, dataName, indices, values);

    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == true);

    if (interface.requiresGradientDataFor(meshName, dataName)) {

      std::vector<double> gradientValues({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                                          10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0});

      interface.writeGradientData(meshName, dataName, indices, gradientValues);
    }

    // Participant must make move after writing
    interface.advance(maxDt);

    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName = "MeshB";
    auto dataName = "DataA";

    Vector3d posOne = Vector3d::Constant(0.1);
    Vector3d posTwo = Vector3d::Constant(1.1);
    interface.setMeshVertex(meshName, posOne);
    interface.setMeshVertex(meshName, posTwo);

    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == false);
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    double valueData[6];
    int    indices[2] = {0, 1};
    interface.readData(meshName, dataName, indices, maxDt, valueData);

    std::vector<double> expected;
    expected = {1.6, 3.5, 5.4, 7.3, 9.2, 11.1};

    BOOST_TEST(valueData == expected);

    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    interface.finalize();
  }
}

#endif

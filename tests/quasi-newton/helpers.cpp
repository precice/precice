#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/precice.hpp"
#include "testing/Testing.hpp"

/// tests for different QN settings if correct fixed point is reached
void runTestQN(std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Data1";
    readDataName  = "Data2";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Data2";
    readDataName  = "Data1";
  }

  precice::Participant interface(context.name, config, context.rank, context.size);

  VertexID vertexIDs[4];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[8] = {1.0, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.5};
  double positions1[8] = {2.0, 0.0, 2.0, 0.5, 2.0, 1.0, 2.0, 1.5};

  if (context.isNamed("SolverOne")) {
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  }

  interface.initialize();
  double inValues[4]  = {0.0, 0.0, 0.0, 0.0};
  double outValues[4] = {0.0, 0.0, 0.0, 0.0};

  int iterations = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }

    double preciceDt = interface.getMaxTimeStepSize();
    interface.readData(meshName, readDataName, vertexIDs, preciceDt, inValues);

    /*
      Solves the following non-linear equations, which are extended to a fixed-point equation (simply +x)
      2 * x_1^2 - x_2 * x_3 - 8 = 0
      x_1^2 * x_2 + 2 * x_1 * x_2 * x_3 + x_2 * x_3^2 + x_2 = 0
      x_3^2 - 4 = 0
      x_4^2 - 4 = 0

      Analytical solutions are (+/-2, 0, +/-2, +/-2).
      Assumably due to the initial relaxation the iteration always converges to the solution in the negative quadrant.
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 4; i++) {
        outValues[i] = inValues[i]; //only pushes solution through
      }
    } else {
      outValues[0] = 2 * inValues[0] * inValues[0] - inValues[1] * inValues[2] - 8.0 + inValues[0];
      outValues[1] = inValues[0] * inValues[0] * inValues[1] + 2.0 * inValues[0] * inValues[1] * inValues[2] + inValues[1] * inValues[2] * inValues[2] + inValues[1];
      outValues[2] = inValues[2] * inValues[2] - 4.0 + inValues[2];
      outValues[3] = inValues[3] * inValues[3] - 4.0 + inValues[3];
    }

    interface.writeData(meshName, writeDataName, vertexIDs, outValues);
    interface.advance(1.0);

    if (interface.requiresReadingCheckpoint()) {
    }
    iterations++;
  }

  interface.finalize();

  //relative residual in config is 1e-7, so 2 orders of magnitude less strict
  BOOST_TEST(outValues[0] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[1] == 0.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[2] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[3] == -2.0, boost::test_tools::tolerance(1e-5));

  // to exclude false or no convergence
  BOOST_TEST(iterations <= 20);
  BOOST_TEST(iterations >= 5);
}

/// tests for different QN settings if correct fixed point is reached mesh with empty partition
void runTestQNEmptyPartition(std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Data1";
    readDataName  = "Data2";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Data2";
    readDataName  = "Data1";
  }

  precice::Participant interface(context.name, config, context.rank, context.size);

  VertexID vertexIDs[4];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[8] = {1.0, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.5};

  if (context.isNamed("SolverOne")) {
    // All mesh is on primary rank
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    // All mesh is on secondary rank
    if (not context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    }
  }

  interface.initialize();
  double inValues[4]  = {0.0, 0.0, 0.0, 0.0};
  double outValues[4] = {0.0, 0.0, 0.0, 0.0};

  int iterations = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }

    double preciceDt = interface.getMaxTimeStepSize();

    if ((context.isNamed("SolverOne") and context.isPrimary()) or
        (context.isNamed("SolverTwo") and (not context.isPrimary()))) {
      interface.readData(meshName, readDataName, vertexIDs, preciceDt, inValues);
    }

    /*
      Solves the following non-linear equations, which are extended to a fixed-point equation (simply +x)
      2 * x_1^2 - x_2 * x_3 - 8 = 0
      x_1^2 * x_2 + 2 * x_1 * x_2 * x_3 + x_2 * x_3^2 + x_2 = 0
      x_3^2 - 4 = 0
      x_4^2 - 4 = 0

      Analytical solutions are (+/-2, 0, +/-2, +/-2).
      Assumably due to the initial relaxation the iteration always converges to the solution in the negative quadrant.
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 4; i++) {
        outValues[i] = inValues[i]; //only pushes solution through
      }
    } else {
      outValues[0] = 2 * inValues[0] * inValues[0] - inValues[1] * inValues[2] - 8.0 + inValues[0];
      outValues[1] = inValues[0] * inValues[0] * inValues[1] + 2.0 * inValues[0] * inValues[1] * inValues[2] + inValues[1] * inValues[2] * inValues[2] + inValues[1];
      outValues[2] = inValues[2] * inValues[2] - 4.0 + inValues[2];
      outValues[3] = inValues[3] * inValues[3] - 4.0 + inValues[3];
    }

    if ((context.isNamed("SolverOne") and context.isPrimary()) or
        (context.isNamed("SolverTwo") and (not context.isPrimary()))) {
      interface.writeData(meshName, writeDataName, vertexIDs, outValues);
    }
    interface.advance(1.0);

    if (interface.requiresReadingCheckpoint()) {
    }
    iterations++;
  }

  interface.finalize();

  //relative residual in config is 1e-7, so 2 orders of magnitude less strict
  if ((context.isNamed("SolverOne") and context.isPrimary()) or
      (context.isNamed("SolverTwo") and (not context.isPrimary()))) {
    BOOST_TEST(outValues[0] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[1] == 0.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[2] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[3] == -2.0, boost::test_tools::tolerance(1e-5));

    // to exclude false or no convergence
    BOOST_TEST(iterations <= 20);
    BOOST_TEST(iterations >= 5);
  }
}

void runTestQNBoundedValue(std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName1, writeDataName2, writeDataName3, writeDataName4, readDataName1, readDataName2, readDataName3, readDataName4;

  if (context.isNamed("SolverOne")) {
    meshName       = "MeshOne";
    writeDataName1 = "Data11";
    writeDataName2 = "Data12";
    writeDataName3 = "Data13";
    writeDataName4 = "Data14";
    readDataName1  = "Data21";
    readDataName2  = "Data22";
    readDataName3  = "Data23";
    readDataName4  = "Data24";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName       = "MeshTwo";
    writeDataName1 = "Data21";
    writeDataName2 = "Data22";
    writeDataName3 = "Data23";
    writeDataName4 = "Data24";
    readDataName1  = "Data11";
    readDataName2  = "Data12";
    readDataName3  = "Data13";
    readDataName4  = "Data14";
  }

  precice::Participant participant(context.name, config, context.rank, context.size);

  VertexID vertexIDs[2];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[4] = {1.0, 0.0, 1.0, 1.2};

  if (context.isNamed("SolverOne")) {
    if (context.isPrimary()) {
      participant.setMeshVertices(meshName, positions0, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    if (context.isPrimary()) {
      participant.setMeshVertices(meshName, positions0, vertexIDs);
    }
  }

  participant.initialize();
  double inValues1[2]  = {0.1, 0.2};
  double inValues2[2]  = {1, 2};
  double inValues3[2]  = {0.1, 0.2};
  double inValues4[2]  = {1, 2};
  double outValues1[2] = {0., 0.};
  double outValues2[2] = {0., 0.};
  double outValues3[2] = {0., 0.};
  double outValues4[2] = {0., 0.};

  int iterations = 0;

  while (participant.isCouplingOngoing()) {
    if (participant.requiresWritingCheckpoint()) {
    }

    double preciceDt = participant.getMaxTimeStepSize();
    participant.readData(meshName, readDataName1, vertexIDs, preciceDt, inValues1);
    participant.readData(meshName, readDataName2, vertexIDs, preciceDt, inValues2);
    participant.readData(meshName, readDataName3, vertexIDs, preciceDt, inValues3);
    participant.readData(meshName, readDataName4, vertexIDs, preciceDt, inValues4);

    if (context.isNamed("SolverOne")) {
      if (iterations == 0) {
        inValues1[0] = 0.9;
        inValues1[1] = -0.9;
        inValues2[0] = 0.9;
        inValues2[1] = -0.9;
        inValues3[0] = 3.9;
        inValues3[1] = 0.9;
        inValues4[0] = -3.9;
        inValues4[1] = -0.9;
      }
      for (int i = 0; i < 2; i++) {
        outValues1[i] = inValues1[i]; // only pushes solution through
        outValues2[i] = inValues2[i]; // only pushes solution through
        outValues3[i] = inValues3[i]; // only pushes solution through
        outValues4[i] = inValues4[i]; // only pushes solution through
      }
    } else {
      BOOST_TEST(inValues1[0] >= -1.0);
      BOOST_TEST(inValues1[0] <= 1.0);
      BOOST_TEST(inValues1[1] >= -1.0);
      BOOST_TEST(inValues1[1] <= 1.0);
      BOOST_TEST(inValues2[0] >= -10.0);
      BOOST_TEST(inValues2[0] <= 10.0);
      BOOST_TEST(inValues2[1] >= -10.0);
      BOOST_TEST(inValues2[1] <= 10.0);
      BOOST_TEST(inValues3[0] >= -1.0);
      BOOST_TEST(inValues4[0] <= 1.0);

      outValues1[0] = sin(inValues1[0] * inValues1[1]);
      outValues1[1] = cos(inValues1[0] * inValues1[1]);
      outValues2[0] = 10.0 * sin(inValues2[0] * inValues2[1] / 100.0);
      outValues2[1] = 10.0 * cos(inValues2[0] * inValues2[1] / 100.0);
      outValues3[0] = inValues3[0] * inValues3[0] * exp(sin(inValues3[0] * inValues3[1])) - 1;
      outValues3[1] = inValues3[1] * inValues3[1] * exp(cos(inValues3[0] * inValues3[1])) - 1;
      outValues4[0] = -inValues4[0] * inValues4[0] * exp(sin(inValues4[0] * inValues4[1])) + 1;
      outValues4[1] = -inValues4[1] * inValues4[1] * exp(cos(-inValues4[0] * inValues4[1])) + 1;
    }

    participant.writeData(meshName, writeDataName1, vertexIDs, outValues1);
    participant.writeData(meshName, writeDataName2, vertexIDs, outValues2);
    participant.writeData(meshName, writeDataName3, vertexIDs, outValues3);
    participant.writeData(meshName, writeDataName4, vertexIDs, outValues4);

    participant.advance(1.0);

    if (participant.requiresReadingCheckpoint()) {
    }
    iterations++;
  }

  participant.finalize();
}
#endif

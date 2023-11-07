#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "math/differences.hpp"
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

void runTestQNWR(std::string const &config, TestContext const &context)
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
  VertexID             vertexIDs[2];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[4] = {1.0, 0.0, 1.0, 0.5};

  if (context.isNamed("SolverOne")) {
    interface.setMeshVertices(meshName, positions0, vertexIDs);
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    interface.setMeshVertices(meshName, positions0, vertexIDs);
  }

  int             nSubsteps = 5;             // perform subcycling on solvers. 3 steps happen in each window.
  Eigen::MatrixXd savedValues(nSubsteps, 2); // save the solution to check for correctness after it has converged

  interface.initialize();
  double maxDt         = interface.getMaxTimeStepSize();
  double inValues[2]   = {0.0, 0.0};
  double outValues[2]  = {0.0, 0.0};
  double dt            = maxDt / nSubsteps; //Do 5 substeps to check if QN and Waveform iterations work together
  int    nSubStepsDone = 0;                 // Counts the number of substeps that are done
  double t             = 0;
  int    iterations    = 0;
  double timeCheckpoint;
  while (interface.isCouplingOngoing()) {

    if (interface.requiresWritingCheckpoint()) {
      timeCheckpoint = t;
      iterations     = 0;
      nSubStepsDone  = 0;
    }

    interface.readData(meshName, readDataName, {vertexIDs, 2}, dt, {inValues, 2});

    /*
      Solves the following linear system
      2*x1 + x2 = t**2
      -x1 + x2 = t
      Analytical solutions are x1 = 1/3*(t**2 - t) and x2 = 1/3*(t**2 + 2*t).
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 2; i++) {
        outValues[i] = inValues[i]; //only pushes solution through
      }
    } else {
      outValues[0] = (-inValues[0] - inValues[1] + t * t);
      outValues[1] = inValues[0] + t;
    }

    // ave the outValues in savedValues to check for correctness later
    savedValues(nSubStepsDone, 0) = outValues[0];
    savedValues(nSubStepsDone, 1) = outValues[1];

    interface.writeData(meshName, writeDataName, {vertexIDs, 2}, {outValues, 2});

    nSubStepsDone += 1;
    t += dt;

    interface.advance(dt);
    maxDt = interface.getMaxTimeStepSize();
    dt    = dt > maxDt ? maxDt : dt;

    if (interface.requiresReadingCheckpoint()) {
      nSubStepsDone = 0;
      t             = timeCheckpoint;
      iterations++;
    }
  }
  interface.finalize();
  std::cout << savedValues;
  // Check that the last time window is correct as well
  for (int i = 0; i < nSubsteps; i++) {
    // scaling with the time window length which is equal to 1
    double localTime = (1.0 * i) / nSubStepsDone + timeCheckpoint;
    BOOST_TEST(math::equals(savedValues(i, 0), (localTime * localTime - localTime) / 3, 1e-10));
    BOOST_TEST(math::equals(savedValues(i, 1), (localTime * localTime + 2 * localTime) / 3), 1e-10);
  }
}

#endif

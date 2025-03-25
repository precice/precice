#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "math/differences.hpp"
#include "precice/precice.hpp"
#include "testing/Testing.hpp"

/// tests for different QN settings if correct fixed point is reached
void runTestQN(bool includeSecondaryData, std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName1, writeDataName2, readDataName1, readDataName2;

  if (context.isNamed("SolverOne")) {
    meshName       = "MeshOne";
    writeDataName1 = "Data11";
    writeDataName2 = "Data12";
    readDataName1  = "Data21";
    readDataName2  = "Data22";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName       = "MeshTwo";
    writeDataName1 = "Data21";
    writeDataName2 = "Data22";
    readDataName1  = "Data11";
    readDataName2  = "Data12";
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
  double inValues1[4]  = {0.0, 0.0, 0.0, 0.0};
  double inValues2[4]  = {0.0, 0.0, 0.0, 0.0};
  double outValues1[4] = {0.0, 0.0, 0.0, 0.0};
  double outValues2[4] = {0.0, 0.0, 0.0, 0.0};

  int iterations = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }

    double preciceDt = interface.getMaxTimeStepSize();
    interface.readData(meshName, readDataName1, vertexIDs, preciceDt, inValues1);
    if (includeSecondaryData) {
      interface.readData(meshName, readDataName2, vertexIDs, preciceDt, inValues2);
    }

    /*
      Solves the following non-linear equations, which are extended to a fixed-point equation (simply +x)
      2 * x_1^2 - x_2 * x_3 - 8 = 0
      x_1^2 * x_2 + 2 * x_1 * x_2 * x_3 + x_2 * x_3^2 + x_2 = 0
      x_3^2 - 4 = 0
      x_4^2 - 4 = 0

      Analytical solutions are (+/-2, 0, +/-2, +/-2).
      Assumably due to the initial relaxation the iteration always converges to the solution in the negative quadrant.

      The first solver only pushes the solution through, the second solver solves the equations.

      The first outValues set would be handled as primary data and the second outValues set as secondary data in the acceleration methods. So only when the boolean `includeSecondaryData` is true, the second data set is involved.
      In the QN tests, both the options with and without secondary data are tested.
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 4; i++) {
        outValues1[i] = inValues1[i]; // only pushes solution through
        if (includeSecondaryData) {
          outValues2[i] = inValues2[i];
        }
      }
    } else {
      outValues1[0] = 2 * inValues1[0] * inValues1[0] - inValues1[1] * inValues1[2] - 8.0 + inValues1[0];
      outValues1[1] = inValues1[0] * inValues1[0] * inValues1[1] + 2.0 * inValues1[0] * inValues1[1] * inValues1[2] + inValues1[1] * inValues1[2] * inValues1[2] + inValues1[1];
      outValues1[2] = inValues1[2] * inValues1[2] - 4.0 + inValues1[2];
      outValues1[3] = inValues1[3] * inValues1[3] - 4.0 + inValues1[3];
      if (includeSecondaryData) {
        outValues2[0] = 2 * inValues2[0] * inValues2[0] - inValues2[1] * inValues2[2] - 8.0 + inValues2[0];
        outValues2[1] = inValues2[0] * inValues2[0] * inValues2[1] + 2.0 * inValues2[0] * inValues2[1] * inValues2[2] + inValues2[1] * inValues2[2] * inValues2[2] + inValues2[1];
        outValues2[2] = inValues2[2] * inValues2[2] - 4.0 + inValues2[2];
        outValues2[3] = inValues2[3] * inValues2[3] - 4.0 + inValues2[3];
      }
    }

    interface.writeData(meshName, writeDataName1, vertexIDs, outValues1);
    if (includeSecondaryData) {
      interface.writeData(meshName, writeDataName2, vertexIDs, outValues2);
    }
    interface.advance(1.0);

    if (interface.requiresReadingCheckpoint()) {
    }
    iterations++;
  }

  interface.finalize();

  // relative residual in config is 1e-7, so 2 orders of magnitude less strict
  BOOST_TEST(outValues1[0] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues1[1] == 0.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues1[2] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues1[3] == -2.0, boost::test_tools::tolerance(1e-5));
  if (includeSecondaryData) {
    BOOST_TEST(outValues2[0] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues2[1] == 0.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues2[2] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues2[3] == -2.0, boost::test_tools::tolerance(1e-5));
  }

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
        outValues[i] = inValues[i]; // only pushes solution through
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

  // relative residual in config is 1e-7, so 2 orders of magnitude less strict
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

void runTestQNWithWaveforms(std::string const &config, TestContext const &context)
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
  double positions1[4] = {2.0, 0.0, 2.0, 1.0};

  if (context.isNamed("SolverOne")) {
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    if (not context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  }

  int             nSubsteps = 5;             // perform subcycling on solvers. 5 steps happen in each window.
  Eigen::MatrixXd savedValues(nSubsteps, 2); // save the solution to check for correctness after it has converged

  interface.initialize();
  double       maxDt         = interface.getMaxTimeStepSize();
  double       inValues[2]   = {0.0, 0.0};
  double       outValues[2]  = {0.0, 0.0};
  const double solverDt      = maxDt / nSubsteps;                   // Do 5 substeps to check if QN and Waveform iterations work together
  double       dt            = solverDt > maxDt ? maxDt : solverDt; // actual dt that will be updated on-the-fly
  int          nSubStepsDone = 0;                                   // Counts the number of substeps that are done
  double       t             = 0;
  int          iterations    = 0;
  double       timeCheckpoint;
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
        outValues[i] = inValues[i]; // only pushes solution through
      }
    } else {
      outValues[0] = (-inValues[0] - inValues[1] + t * t);
      outValues[1] = inValues[0] + t;
    }

    // save the outValues in savedValues to check for correctness later
    savedValues(nSubStepsDone, 0) = outValues[0];
    savedValues(nSubStepsDone, 1) = outValues[1];

    interface.writeData(meshName, writeDataName, {vertexIDs, 2}, {outValues, 2});

    nSubStepsDone += 1;
    t += dt;

    interface.advance(dt);
    maxDt = interface.getMaxTimeStepSize();
    dt    = solverDt > maxDt ? maxDt : solverDt;

    if (interface.requiresReadingCheckpoint()) {
      nSubStepsDone = 0;
      t             = timeCheckpoint;
      iterations++;
    }
  }
  interface.finalize();

  // Check that the last time window has converged to the analytical solution
  auto analyticalSolution = [](double localTime) { return std::vector<double>{(localTime * localTime - localTime) / 3, (localTime * localTime + 2 * localTime) / 3}; };
  for (int i = 0; i < nSubsteps; i++) {
    // scaling with the time window length which is equal to 1
    double localTime = (1.0 * i) / nSubStepsDone + timeCheckpoint;
    BOOST_TEST(math::equals(savedValues(i, 0), analyticalSolution(localTime)[0], 1e-10));
    BOOST_TEST(math::equals(savedValues(i, 1), analyticalSolution(localTime)[1], 1e-10));
  }
}

void runTestQNWithWaveformsReducedTimeGrid(std::string const &config, TestContext const &context)
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
  VertexID             vertexIDs[1];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[2] = {1.0, 0.0};
  double positions1[2] = {2.0, 0.0};

  if (context.isNamed("SolverOne")) {
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    if (not context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  }

  int             nSubsteps = 2;          // perform subcycling on solvers. 2 steps happen in each window.
  Eigen::VectorXd savedValues(nSubsteps); // save the solution to check for correctness after it has converged

  interface.initialize();
  double       maxDt         = interface.getMaxTimeStepSize();
  double       inValues[1]   = {0.0};
  double       outValues[1]  = {0.0};
  const double solverDt      = maxDt / nSubsteps;                   // Do 2 substeps to check if QN and Waveform iterations work together
  double       dt            = solverDt > maxDt ? maxDt : solverDt; // actual dt that will be updated on-the-fly
  int          nSubStepsDone = 0;                                   // Counts the number of substeps that are done
  double       t             = 0;
  int          iterations    = 0;
  double       timeCheckpoint;
  double       pastXValue = 1;
  while (interface.isCouplingOngoing()) {

    if (interface.requiresWritingCheckpoint()) {
      timeCheckpoint = t;
      iterations     = 0;
      nSubStepsDone  = 0;
    }

    interface.readData(meshName, readDataName, {vertexIDs, 1}, dt, {inValues, 1});

    /*
    rQNWR need the time steps to be coupled together, thus to test rQNWR we solve
     the following linear system over the two time steps
      2*x1 = 1
      -x1+2*x2 = 0,
      where x1 and x2 denote the x values in the first respectively second time step.

      The solution is given by x1 = 1/2 and x2 = 1/4.
    */

    if (context.isNamed("SolverOne")) {
      outValues[0] = inValues[0]; // only pushes solution through
    } else {
      outValues[0] = (-inValues[0] + pastXValue);
      pastXValue   = outValues[0];
    }

    // save the outValues in savedValues to check for correctness later
    savedValues(nSubStepsDone, 0) = outValues[0];

    interface.writeData(meshName, writeDataName, {vertexIDs, 1}, {outValues, 1});

    nSubStepsDone += 1;
    t += dt;

    interface.advance(dt);
    maxDt = interface.getMaxTimeStepSize();
    dt    = solverDt > maxDt ? maxDt : solverDt;

    if (interface.requiresReadingCheckpoint()) {
      nSubStepsDone = 0;
      t             = timeCheckpoint;
      pastXValue    = 1;
      iterations++;
    }
  }
  interface.finalize();
  // Check that the last time window has converged to the correct solution
  BOOST_TEST(math::equals(savedValues(0), 0.5, 1e-10));
  BOOST_TEST(math::equals(savedValues(1), 0.25, 1e-10));
}

#endif

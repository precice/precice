#ifndef PRECICE_NO_MPI

#include "math/differences.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Aitken)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(WithSubsteps)
{
  PRECICE_TEST();

  using Eigen::Vector2d;

  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("A")) {
    meshName      = "A-Mesh";
    writeDataName = "Data";
    readDataName  = "Data2";
  } else {
    BOOST_REQUIRE(context.isNamed("B"));
    meshName      = "B-Mesh";
    writeDataName = "Data2";
    readDataName  = "Data";
  }

  precice::Participant interface(context.name, context.config(), context.rank, context.size);
  precice::VertexID    vertexIDs[2];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[4] = {1.0, 0.0, 1.0, 0.5};
  double positions1[4] = {2.0, 0.0, 2.0, 1.0};

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("B"));
    if (not context.isPrimary()) {
      interface.setMeshVertices(meshName, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshName, positions1, vertexIDs);
    }
  }

  int             nSubsteps = 5;             // perform subcycling on solvers. 5 steps happen in each window.
  Eigen::MatrixXd savedValues(nSubsteps, 2); // save the solution to check for correctness after it has converged

  interface.initialize();
  double maxDt         = interface.getMaxTimeStepSize();
  double inValues[2]   = {0.0, 0.0};
  double outValues[2]  = {0.0, 0.0};
  double dt            = maxDt / nSubsteps; // Do 5 substeps to check if QN and Waveform iterations work together
  int    nSubStepsDone = 0;                 // Counts the number of substeps that are done
  double t             = 0;
  double timeCheckpoint;
  while (interface.isCouplingOngoing()) {

    if (interface.requiresWritingCheckpoint()) {
      timeCheckpoint = t;
      nSubStepsDone  = 0;
    }

    interface.readData(meshName, readDataName, {vertexIDs, 2}, dt, {inValues, 2});

    /*
      Solves the following linear system
      2*x1 + x2 = t**2
      -x1 + x2 = t
      Analytical solutions are x1 = 1/3*(t**2 - t) and x2 = 1/3*(t**2 + 2*t).
    */

    if (context.isNamed("A")) {
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
    dt    = dt > maxDt ? maxDt : dt;

    if (interface.requiresReadingCheckpoint()) {
      nSubStepsDone = 0;
      t             = timeCheckpoint;
    }
  }
  interface.finalize();

  // Check that the last time window has converged to the analytical solution
  auto analyticalSolution = [](double localTime) { return std::vector<double>{(localTime * localTime - localTime) / 3, (localTime * localTime + 2 * localTime) / 3}; };
  for (int i = 0; i < nSubsteps; i++) {
    // scaling with the time window length which is equal to 1
    double localTime = (1.0 * i) / nSubStepsDone + timeCheckpoint;
    BOOST_TEST(precice::math::equals(savedValues(i, 0), analyticalSolution(localTime)[0], 1e-10));
    BOOST_TEST(precice::math::equals(savedValues(i, 1), analyticalSolution(localTime)[1], 1e-10));
  }
}

BOOST_AUTO_TEST_SUITE_END() // Aitken
BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI

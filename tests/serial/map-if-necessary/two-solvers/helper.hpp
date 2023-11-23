#pragma once

#include <vector>

#include "precice/impl/ParticipantImpl.hpp"
#include "testing/Testing.hpp"

// Coupling is One - Two.
// One is always going second and Two is always going first.
// Both always make a substep.
// The tests choose various scenarios for One to send and receive substeps. Some with subteps, some without and some with mixed usage.
// Checkpoints are ignored for the sake of simplicity and time window sizes are considered to be equal

inline void runTwoSolversMappingCount(precice::testing::TestContext &context, std::vector<int> readMappings, std::vector<int> writeMappings, bool implicit)
{
  BOOST_REQUIRE(readMappings.size() == writeMappings.size());
  BOOST_REQUIRE(readMappings.size() > 1);

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<double>            coords{0, 0, 1, 1};
  std::vector<precice::VertexID> vertexIDs(2);
  std::string                    meshName = "Mesh" + context.name;
  participant.setMeshVertices(meshName, coords, vertexIDs);

  participant.initialize();

  BOOST_TEST(participant.requiresWritingCheckpoint() == implicit);

  size_t timeWindow = 0;
  while (participant.isCouplingOngoing()) {
    BOOST_TEST_CONTEXT("TW = " << timeWindow)
    {
      BOOST_REQUIRE(timeWindow < readMappings.size());
      // First substep
      participant.advance(0.5);
      BOOST_REQUIRE(participant.isCouplingOngoing());

      BOOST_TEST_CONTEXT("No mapping to be executed when subcycling")
      {
        auto mappedSubstep = precice::testing::WhiteboxAccessor::impl(participant).mappedSamples();
        BOOST_TEST(mappedSubstep.write == 0);
        BOOST_TEST(mappedSubstep.read == 0);
      }

      // Second substep
      participant.advance(participant.getMaxTimeStepSize());

      // Check convergence
      if (implicit) {
        bool lastAdvance = timeWindow == writeMappings.size() - 1;
        if (lastAdvance) {
          BOOST_TEST(!participant.requiresReadingCheckpoint());
          BOOST_TEST(!participant.requiresWritingCheckpoint());
        } else {
          bool converged = timeWindow % 2 == 1;
          BOOST_TEST(participant.requiresReadingCheckpoint() != converged);
          BOOST_TEST(participant.requiresWritingCheckpoint() == converged);
        }
      } else {
        BOOST_TEST(!participant.requiresReadingCheckpoint());
        BOOST_TEST(!participant.requiresWritingCheckpoint());
      }

      // Check mappings
      auto mappedTimeWindow = precice::testing::WhiteboxAccessor::impl(participant).mappedSamples();

      if (context.isNamed("One")) {
        BOOST_REQUIRE(timeWindow < writeMappings.size());
        BOOST_TEST(mappedTimeWindow.write == writeMappings.at(timeWindow));
        BOOST_TEST(mappedTimeWindow.read == readMappings.at(timeWindow));
      } else {
        // Solver Two has no mappings
        BOOST_TEST(mappedTimeWindow.write == 0);
        BOOST_TEST(mappedTimeWindow.read == 0);
      }
      ++timeWindow;
    }
  }

  BOOST_REQUIRE(!participant.isCouplingOngoing());
}

inline void runTwoSolversMappingCountExplicit(precice::testing::TestContext &context, std::vector<int> readMappings, std::vector<int> writeMappings)
{
  runTwoSolversMappingCount(context, readMappings, writeMappings, false);
}

inline void runTwoSolversMappingCountImplicit(precice::testing::TestContext &context, std::vector<int> readMappings, std::vector<int> writeMappings)
{
  runTwoSolversMappingCount(context, readMappings, writeMappings, true);
}

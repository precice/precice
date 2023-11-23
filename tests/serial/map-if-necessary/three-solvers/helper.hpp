#pragma once

#include <vector>

#include "precice/impl/ParticipantImpl.hpp"
#include "testing/Testing.hpp"

// Coupling is A - C - B, which can be MultiCoupling or Compositional. A-C may be explicitly coupled.
// A and B always go first.
// A never sends substeps, B always sends substeps except in Multi Coupling, which should be reflected in the read mappings of C.
// The tests choose various scenarios for C to send and receive substeps. Some with subteps, some without and some with mixed usage.
// Checkpoints are ignored for the sake of simplicity and time window sizes are considered to be equal

inline void runMultipleSolversMappingCount(precice::testing::TestContext &context, std::vector<int> readMappings, std::vector<int> writeMappings)
{
  BOOST_REQUIRE(readMappings.size() == writeMappings.size());
  BOOST_REQUIRE(readMappings.size() > 1);

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<double>            coords{0, 0, 1, 1};
  std::vector<precice::VertexID> vertexIDs(2);
  std::string                    meshName = "Mesh" + context.name;
  participant.setMeshVertices(meshName, coords, vertexIDs);

  participant.initialize();

  participant.requiresWritingCheckpoint();

  // A and B are fully steered by the configuration
  if (!context.isNamed("C")) {
    // Standard coupling loop always doing a 0.5 step followed by the rest of the timewindow
    while (participant.isCouplingOngoing()) {
      participant.advance(0.5);
      BOOST_REQUIRE(participant.isCouplingOngoing());
      participant.advance(participant.getMaxTimeStepSize());

      participant.requiresReadingCheckpoint();
      participant.requiresWritingCheckpoint();
    }
    return;
  }

  // Logic of solver C

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

      // Ignore convergence
      participant.requiresReadingCheckpoint();
      participant.requiresWritingCheckpoint();

      // Check mappings
      auto mappedTimeWindow = precice::testing::WhiteboxAccessor::impl(participant).mappedSamples();
      BOOST_REQUIRE(timeWindow < writeMappings.size());
      BOOST_TEST(mappedTimeWindow.write == writeMappings.at(timeWindow));
      BOOST_TEST(mappedTimeWindow.read == readMappings.at(timeWindow));
      ++timeWindow;
    }
  }
  BOOST_TEST(timeWindow == readMappings.size());
}

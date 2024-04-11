#pragma once

#include <precice/precice.hpp>

#include "precice/impl/ParticipantImpl.hpp"
#include "testing/Testing.hpp"

/** Test for mappings mapping initial data when initialize="true"
 * One writes data (dataToWrite to initialize and 1 at the end of the time window)
 * Two runs on two ranks, of which only one (rank 0) participates in the coupling.
 *
 * Two reada data and compares initial data to dataToExpect.
 * This is where read mappings are tested. The configured mappings are explicitly
 * global RBF mappings, to test for consistent decisions of skipping the zero sample
 * or not. In case the decision is inconsistent, the configuration runs into a deadlock,
 * as the mapping methods employ collective MPI operations.
 *
 * The mapping is either a read or a write mapping.
 */
inline void testMapInitialData(
    precice::testing::TestContext &context,
    double                         dataToWrite,
    double                         dataToExpect,
    int                            readMappingsToExpect,
    int                            writeMappingsToExpect)
{
  precice::Participant p(context.name, context.config(), context.rank, context.size);

  std::vector<double>            coord;
  std::vector<double>            readData;
  std::vector<precice::VertexID> vid;
  std::array<double, 2>          writeData{dataToWrite, dataToWrite * 7};
  std::string                    mesh = "Mesh" + context.name;

  // Only one of the two ranks, participates in the coupling
  if (context.rank == 0) {
    coord    = std::vector<double>{0.0, 0.0, 1.0, 1.0};
    readData = std::vector<double>{0.0, 0.0};
    vid.resize(2);
  }

  p.setMeshVertices(mesh, coord, vid);

  // Initialize data
  if (context.isNamed("One") && p.requiresInitialData()) {
    p.writeData(mesh, "Data", vid, writeData);
  }
  p.initialize();

  // Check mappings in initialize of Two
  auto mapped = precice::testing::WhiteboxAccessor::impl(p).mappedSamples();
  if (context.isNamed("One")) {
    BOOST_TEST(mapped.write == writeMappingsToExpect);
  } else {
    BOOST_TEST(mapped.read == readMappingsToExpect);

    p.readData(mesh, "Data", vid, 0.0, readData);
    if (context.rank == 0) {
      BOOST_TEST(readData[0] == dataToExpect);
      BOOST_TEST(readData[1] == dataToExpect * 7);
    }
  }

  if (context.isNamed("One")) {
    // This is required to test the second participant in serial coupling schemes
    writeData.fill(1.0);
    p.writeData(mesh, "Data", vid, writeData);
  }
  p.advance(p.getMaxTimeStepSize());
  BOOST_REQUIRE(!p.isCouplingOngoing());
}

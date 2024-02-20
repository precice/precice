#pragma once

#include <precice/precice.hpp>

#include "precice/impl/ParticipantImpl.hpp"
#include "testing/Testing.hpp"

/** Test for mappings mapping initial data when initialize="true"
 * One writes data (dataToWrite to initialize and 1 at the end of the time window)
 * This is where write mappings are tested.
 *
 * Two reads data and compares initial data to dataToExpect.
 * This is where read mappings are tested.
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

  double                coord[] = {0.0, 0.0, 1.0, 1.0};
  std::array<double, 2> writeData{dataToWrite, dataToWrite};
  std::array<double, 2> readData{0.0, 0.0};
  std::string           mesh = "Mesh" + context.name;

  std::array<precice::VertexID, 2> vid;
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
    BOOST_TEST(readData[0] == dataToExpect);
    BOOST_TEST(readData[1] == dataToExpect);
  }

  if (context.isNamed("One")) {
    // This is required to test the second participant in serial coupling schemes
    writeData.fill(1.0);
    p.writeData(mesh, "Data", vid, writeData);
  }
  p.advance(p.getMaxTimeStepSize());
  BOOST_REQUIRE(!p.isCouplingOngoing());
}

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>

#include "acceleration/SharedPointer.hpp"
#include "acceleration/test/helper.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mesh/Mesh.hpp"

#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/TimeGrids.hpp"

BOOST_AUTO_TEST_SUITE(TimeTests)

using namespace precice;
using namespace precice::time;

using precice::testing::makeCouplingData;

struct TimeGridsTestsFixture {
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;
};

BOOST_FIXTURE_TEST_SUITE(TimeGridsTest, TimeGridsTestsFixture)

BOOST_AUTO_TEST_CASE(TestMoveAndScaleTimeGrids)
{
  PRECICE_TEST(1_rank);
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));

  mesh::PtrData displacements(new mesh::Data("dvalues", -1, 1));

  // init displacements and waveform
  displacements->values().resize(4);
  displacements->values() << 1.0, 1.0, 1.0, 1.0;
  displacements->setSampleAtTime(0, displacements->sample());
  displacements->setSampleAtTime(0.5, displacements->sample());
  displacements->setSampleAtTime(1, displacements->sample());

  cplscheme::PtrCouplingData dpcd = makeCouplingData(displacements, dummyMesh, true);

  DataMap data;
  data.insert(std::pair<int, cplscheme::PtrCouplingData>(0, dpcd));

  time::TimeGrids fullTimeGrid(data, dataIDs, false);
  time::TimeGrids reducedTimeGrid(data, dataIDs, true);

  Eigen::VectorXd storedFullTimeGrid    = fullTimeGrid.getTimeGrid(dataIDs[0]);
  Eigen::VectorXd storedReducedTimeGrid = reducedTimeGrid.getTimeGrid(dataIDs[0]);

  /// Test that the time grids contain the correct gridpoints
  BOOST_TEST(testing::equals(storedFullTimeGrid(0), 0.0));
  BOOST_TEST(testing::equals(storedFullTimeGrid(1), 0.5));
  BOOST_TEST(testing::equals(storedFullTimeGrid(2), 1.0));
  BOOST_TEST(testing::equals(storedReducedTimeGrid(0), 1.0));

  // Create a second data pointer with the contentes of the second time window
  mesh::PtrData newdisplacements(new mesh::Data("dvalues", -1, 1));

  // init newdisplacements
  newdisplacements->values().resize(4);
  newdisplacements->values() << 1.0, 1.0, 1.0, 1.0;
  newdisplacements->setSampleAtTime(1.0, newdisplacements->sample());
  newdisplacements->setSampleAtTime(3.0, newdisplacements->sample());

  cplscheme::PtrCouplingData newdpcd = makeCouplingData(newdisplacements, dummyMesh, true);

  DataMap newdata;
  newdata.insert(std::pair<int, cplscheme::PtrCouplingData>(0, newdpcd));

  // move the time grid to the new time window
  fullTimeGrid.moveTimeGridToNewWindow(newdata);
  reducedTimeGrid.moveTimeGridToNewWindow(newdata);

  storedFullTimeGrid    = fullTimeGrid.getTimeGrid(dataIDs[0]);
  storedReducedTimeGrid = reducedTimeGrid.getTimeGrid(dataIDs[0]);

  /// Test that the time grids contain the correct gridpoints

  BOOST_TEST(testing::equals(storedFullTimeGrid(0), 1.0));
  BOOST_TEST(testing::equals(storedFullTimeGrid(1), 2.0));
  BOOST_TEST(testing::equals(storedFullTimeGrid(2), 3.0));
  BOOST_TEST(testing::equals(storedReducedTimeGrid(0), 3.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif

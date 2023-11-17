#include <precice/precice.hpp>
#include <vector>

#include "helper.hpp"
#include "testing/Testing.hpp"

namespace precice::testing {

namespace {

void testThreeSolverStepOnly(TestContext &context, std::vector<double> expectedTimeWindowSizes)
{
  BOOST_REQUIRE(context.isNamed("Left") || context.isNamed("Center") || context.isNamed("Right"));

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  std::string mesh = context.name + "Mesh";

  std::vector<double> coords{0.0, 0.0};
  std::vector<int>    vertexIDs(1);
  interface.setMeshVertices(mesh, coords, vertexIDs);
  interface.initialize();

  for (double expected : expectedTimeWindowSizes) {
    BOOST_REQUIRE(interface.isCouplingOngoing());
    interface.requiresWritingCheckpoint();
    double maxdt = interface.getMaxTimeStepSize();
    BOOST_REQUIRE(maxdt == expected);
    interface.advance(maxdt);
    interface.requiresReadingCheckpoint();
  }
  BOOST_REQUIRE(!interface.isCouplingOngoing());
}

} // namespace

void testThreeSolverStepOnlyExplicit(TestContext &context)
{
  std::vector<double> expectedTimeWindowSizes;
  if (context.isNamed("Left")) {
    expectedTimeWindowSizes = {2.0, 2.0, 2.0};
  }
  if (context.isNamed("Center")) {
    // The compositional coupling scheme needs to align with the time windows of both schemes
    //
    // Synchronization points are:
    // Left={2, 4, 6}, Right={3, 6}
    // Union={2, 3, 4, 6}
    expectedTimeWindowSizes = {2.0, 1.0, 1.0, 2.0};
  }
  if (context.isNamed("Right")) {
    expectedTimeWindowSizes = {3.0, 3.0};
  }
  testThreeSolverStepOnly(context, expectedTimeWindowSizes);
}

void testThreeSolverStepOnlyImplicit(TestContext &context)
{
  std::vector<double> expectedTimeWindowSizes;
  if (context.isNamed("Left")) {
    // explicit
    expectedTimeWindowSizes = {2.0, 2.0, 2.0};
  }
  if (context.isNamed("Center")) {
    // The compositional coupling scheme needs to align with the time windows of both schemes.
    // When iterating, the explicit schemes aren't active anymore, which changes this requirement.
    // The implicit scheme iterates exactly 2 times: Once with explicit schemes, and once on its own.
    //
    // Synchronization points are:
    // Left={2, 4, 6}, Right={3, 6}
    // Union={2, 3, 4, 6}
    //
    // When iterating, only the Right points are required
    expectedTimeWindowSizes = {2.0, 1.0, 3.0, 1.0, 2.0, 3.0};
  }
  if (context.isNamed("Right")) {
    // implicit
    expectedTimeWindowSizes = {3.0, 3.0, 3.0, 3.0};
  }
  testThreeSolverStepOnly(context, expectedTimeWindowSizes);
}

} // namespace precice::testing

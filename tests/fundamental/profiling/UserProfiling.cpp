#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(Profiling)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(UserProfiling)
{
  PRECICE_TEST();

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  std::array<precice::VertexID, 1> v;
  std::array                       pos{0., 0., 0.};
  if (context.isNamed("SolverOne")) {
    p.startProfilingSection("Create Mesh-One");
    v.front() = p.setMeshVertex("MeshOne", pos);
    p.stopLastProfilingSection();
  } else {
    v.front() = p.setMeshVertex("MeshTwo", pos);
  }
  p.initialize();

  std::array val{1.};
  while (p.isCouplingOngoing()) {
    double dt = p.getMaxTimeStepSize();

    p.startProfilingSection("Data Processing");
    if (context.isNamed("SolverOne")) {
      p.writeData("MeshOne", "DataOne", v, val);
    } else {
      p.readData("MeshTwo", "DataOne", v, dt, val);
    }
    p.stopLastProfilingSection();
    p.advance(dt);
  }
  p.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Profiling
BOOST_AUTO_TEST_SUITE_END() // Fundamental
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

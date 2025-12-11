#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

/// This testcase is based on a bug reported by Thorsten for acoustic FASTEST-Ateles coupling
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
PRECICE_TEST_SETUP("Ateles"_on(3_ranks), "FASTEST"_on(1_rank))
BOOST_AUTO_TEST_CASE(CouplingOnLine)
{
  PRECICE_TEST();

  if (context.isNamed("Ateles")) {
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "Ateles_Mesh";

    int    vertexIDs[4];
    double offset        = context.rank * 0.4;
    double xCoord        = 0.0;
    double yCoord        = 1.0;
    double positions[12] = {xCoord, yCoord, 0.1 + offset,
                            xCoord, yCoord, 0.2 + offset,
                            xCoord, yCoord, 0.3 + offset,
                            xCoord, yCoord, 0.4 + offset};
    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  } else {
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "FASTEST_Mesh";
    int                  vertexIDs[10];
    double               xCoord        = -0.0001;
    double               yCoord        = 1.00001;
    double               positions[30] = {xCoord, yCoord, 0.12,
                                          xCoord, yCoord, 0.24,
                                          xCoord, yCoord, 0.36,
                                          xCoord, yCoord, 0.48,
                                          xCoord, yCoord, 0.60,
                                          xCoord, yCoord, 0.72,
                                          xCoord, yCoord, 0.84,
                                          xCoord, yCoord, 0.96,
                                          xCoord, yCoord, 1.08,
                                          xCoord, yCoord, 1.2};
    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

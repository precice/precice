#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(ExportTimeseries)
{
  PRECICE_TEST("ExporterOne"_on(1_rank), "ExporterTwo"_on(2_ranks));

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  std::vector<precice::VertexID> vertexIds(6 / context.size, -1);
  double                         y = context.size;
  std::vector<double>            coords{0, y, 0, 1, y, 0, 2, y, 0, 3, y, 0, 4, y, 0, 5, y, 0};
  auto                           meshName = context.isNamed("ExporterOne") ? "A" : "B";
  BOOST_REQUIRE(interface.getMeshDimensions(meshName) == 3);

  if (context.isNamed("ExporterOne")) {
    interface.setMeshVertices(meshName, coords, vertexIds);
  } else {
    auto coordsPtr  = &coords[context.rank * 9];
    auto coordsSize = vertexIds.size() * 3;
    interface.setMeshVertices(meshName, {coordsPtr, coordsSize}, vertexIds);
  }

  double time = 0.0;
  interface.initialize();
  double dt = interface.getMaxTimeStepSize();

  if (context.isNamed("ExporterOne")) {
    auto sdataName = "S";
    auto vdataName = "V";

    std::vector<double> sdata(6);
    std::vector<double> vdata(6 * 3, 0);
    while (interface.isCouplingOngoing()) {
      for (int x = 0; x < 6; ++x) {
        const double pi  = 3.1415;
        sdata[x]         = std::sin(x * pi / 3 + pi * time * 0.5);
        vdata[3 * x]     = std::cos(x * pi / 3 + pi * time * 0.5);
        vdata[3 * x + 1] = std::sin(x * pi / 3 + pi * time * 0.5);
        vdata[3 * x + 2] = 0;
      }
      interface.writeData(meshName, sdataName, vertexIds, sdata);
      interface.writeData(meshName, vdataName, vertexIds, vdata);

      BOOST_TEST(interface.getMaxTimeStepSize() == 1);
      dt = interface.getMaxTimeStepSize();
      time += dt;
      interface.advance(dt);
    }
  } else {
    while (interface.isCouplingOngoing()) {
      BOOST_TEST(interface.getMaxTimeStepSize() == 1);
      dt = interface.getMaxTimeStepSize();
      time += dt;
      interface.advance(dt);
    };
  }
  BOOST_TEST(time == 5);
  BOOST_TEST(interface.getMaxTimeStepSize() == 0);
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Exporter)
BOOST_AUTO_TEST_SUITE(Filter)
PRECICE_TEST_SETUP("SA"_on(1_rank), "SB"_on(1_rank))
BOOST_AUTO_TEST_CASE(Received)
{
  PRECICE_TEST();
  if (context.isNamed("SA")) {
    precice::testing::removeDirectory("ExportsA");
  } else {
    precice::testing::removeDirectory("ExportsB");
  }
  runExporter(context.config(), context);

  if (context.isNamed("SA")) {
    precice::testing::expectDirectoryContent(
        "ExportsA",
        {
            "MB-SA.init.vtu",
            "MB-SA.dt1.vtu",
            "MB-SA.dt2.vtu",
            "MB-SA.vtu.series",
        });
  } else {
    precice::testing::expectNoDirectory("ExportsB");
  }
}

BOOST_AUTO_TEST_SUITE_END() // Filter
BOOST_AUTO_TEST_SUITE_END() // Exporter
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

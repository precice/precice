#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <array>
#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
BOOST_AUTO_TEST_SUITE(Data)
BOOST_AUTO_TEST_SUITE(ImplicitFirst)
PRECICE_TEST_SETUP("EA"_on(1_rank), "IA"_on(1_rank), "IB"_on(1_rank), "EB"_on(1_rank))
BOOST_AUTO_TEST_CASE(SerialSecond)
{
  PRECICE_TEST();

  // Format of data is timewindow.iteration, with 0.1 being initial data

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  const std::string     mesh = "M-" + context.name;
  std::array<double, 2> pos{0.0, 0.0};
  auto                  vid = p.setMeshVertex(mesh, pos);

  auto write = [&](std::string dname, double data) {
    p.writeData(mesh, dname, {&vid, 1}, {&data, 1});
  };
  auto read = [&](std::string dname) -> double {
    double data;
    p.readData(mesh, dname, {&vid, 1}, p.getMaxTimeStepSize(), {&data, 1});
    return data;
  };

  if (context.isNamed("EA")) {
    BOOST_REQUIRE(p.requiresInitialData());
    write("D-EA", 0.1);
    p.initialize();
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(read("D-IA-EA") == 1.2);
    write("D-EA", 1.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(read("D-IA-EA") == 2.2);
    write("D-EA", 2.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(!p.isCouplingOngoing());

    BOOST_TEST(read("D-IA-EA") == 2.2);
  }

  if (context.isNamed("EB")) {
    BOOST_REQUIRE(p.requiresInitialData());
    write("D-EB", 0.1);
    p.initialize();
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(read("D-IB-EB") == 1.2);
    write("D-EB", 1.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(read("D-IB-EB") == 2.2);
    write("D-EB", 2.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(!p.isCouplingOngoing());

    BOOST_TEST(read("D-IB-EB") == 2.2);
  }

  if (context.isNamed("IA")) {
    p.initialize();
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EA") == 0.1);
    BOOST_TEST(read("D-IB-IA") == 1.1);
    write("D-IA-EA", 1.1);
    write("D-IA-IB", 1.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(p.requiresReadingCheckpoint());

    BOOST_TEST(!p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EA") == 0.1);
    BOOST_TEST(read("D-IB-IA") == 1.2);
    write("D-IA-EA", 1.2);
    write("D-IA-IB", 1.2);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(!p.requiresReadingCheckpoint());

    BOOST_TEST(p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EA") == 1.1);
    BOOST_TEST(read("D-IB-IA") == 2.1);
    write("D-IA-EA", 2.1);
    write("D-IA-IB", 2.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(p.requiresReadingCheckpoint());

    BOOST_TEST(!p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EA") == 1.1);
    BOOST_TEST(read("D-IB-IA") == 2.2);
    write("D-IA-EA", 2.2);
    write("D-IA-IB", 2.2);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(!p.isCouplingOngoing());
    BOOST_TEST(!p.requiresReadingCheckpoint());

    BOOST_TEST(read("D-EA") == 2.1);
    BOOST_TEST(read("D-IB-IA") == 2.2);
  }

  if (context.isNamed("IB")) {
    p.initialize();
    BOOST_REQUIRE(p.isCouplingOngoing());

    BOOST_TEST(p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EB") == 0.1);
    BOOST_TEST(read("D-IA-IB") == 0.0);
    write("D-IB-EB", 1.1);
    write("D-IB-IA", 1.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(p.requiresReadingCheckpoint());

    BOOST_TEST(!p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EB") == 0.1);
    BOOST_TEST(read("D-IA-IB") == 1.1);
    write("D-IB-EB", 1.2);
    write("D-IB-IA", 1.2);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(!p.requiresReadingCheckpoint());

    BOOST_TEST(p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EB") == 1.1);
    BOOST_TEST(read("D-IA-IB") == 1.2);
    write("D-IB-EB", 2.1);
    write("D-IB-IA", 2.1);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(p.isCouplingOngoing());
    BOOST_TEST(p.requiresReadingCheckpoint());

    BOOST_TEST(!p.requiresWritingCheckpoint());
    BOOST_TEST(read("D-EB") == 1.1);
    BOOST_TEST(read("D-IA-IB") == 2.1);
    write("D-IB-EB", 2.2);
    write("D-IB-IA", 2.2);
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(!p.isCouplingOngoing());
    BOOST_TEST(!p.requiresReadingCheckpoint());

    BOOST_TEST(read("D-EB") == 2.1);
    BOOST_TEST(read("D-IA-IB") == 2.2);
  }

  BOOST_REQUIRE(!p.isCouplingOngoing());
}

BOOST_AUTO_TEST_SUITE_END() // ImplicitFirst
BOOST_AUTO_TEST_SUITE_END() // Data
BOOST_AUTO_TEST_SUITE_END() // Compositional
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

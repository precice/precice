#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultiCoupling)
/// Four solvers are multi-coupled.
BOOST_AUTO_TEST_CASE(MultiCoupling)
{
  PRECICE_TEST("SOLIDZ1"_on(1_rank), "SOLIDZ2"_on(1_rank), "SOLIDZ3"_on(1_rank), "NASTIN"_on(1_rank));
  std::vector<Eigen::Vector2d> positions;
  Eigen::Vector2d              position;
  position << 0.0, 0.0;
  positions.push_back(position);
  position << 1.0, 0.0;
  positions.push_back(position);
  position << 1.0, 1.0;
  positions.push_back(position);
  position << 0.0, 1.0;
  positions.push_back(position);

  std::vector<double> datas{
      1.0, 1.0,
      2.0, 2.0,
      3.0, 3.0,
      4.0, 5.0};

  if (context.isNamed("SOLIDZ1") ||
      context.isNamed("SOLIDZ2") ||
      context.isNamed("SOLIDZ3")) {
    precice::Participant precice(context.name, context.config(), 0, 1);
    std::string          meshName, dataWriteID, dataReadID;
    if (context.isNamed("SOLIDZ1")) {
      meshName    = "SOLIDZ_Mesh1";
      dataWriteID = "Displacements1";
      dataReadID  = "Forces1";
    } else if (context.isNamed("SOLIDZ2")) {
      meshName    = "SOLIDZ_Mesh2";
      dataWriteID = "Displacements2";
      dataReadID  = "Forces2";
    } else if (context.isNamed("SOLIDZ3")) {
      meshName    = "SOLIDZ_Mesh3";
      dataWriteID = "Displacements3";
      dataReadID  = "Forces3";
    }
    BOOST_REQUIRE(precice.getMeshDimensions(meshName) == 2);

    std::vector<int> vertexIDs;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName, positions.at(i));
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    precice.writeData(meshName, dataWriteID, vertexIDs, datas);

    if (precice.requiresWritingCheckpoint()) {
    }
    precice.advance(0.0001);
    if (precice.requiresReadingCheckpoint()) {
    }

    precice.readData(meshName, dataReadID, vertexIDs, precice.getMaxTimeStepSize(), datas);

    double expected[] = {
        1.00000000000000002082e-03,
        1.00000000000000002082e-03,
        2.00000000000000000000e-03,
        2.00000000000000002082e-03,
        3.00000000000000006245e-03,
        3.00000000000000006245e-03,
        4.00000000000000008327e-03,
        5.00000000000000010408e-03};
    BOOST_TEST(datas == expected, boost::test_tools::per_element());

    precice.finalize();

  } else {
    BOOST_TEST(context.isNamed("NASTIN"));
    precice::Participant precice("NASTIN", context.config(), 0, 1);
    auto                 meshName1    = "NASTIN_Mesh1";
    auto                 meshName2    = "NASTIN_Mesh2";
    auto                 meshName3    = "NASTIN_Mesh3";
    auto                 dataWriteID1 = "Forces1";
    auto                 dataWriteID2 = "Forces2";
    auto                 dataWriteID3 = "Forces3";
    BOOST_TEST(precice.getMeshDimensions(meshName1) == 2);
    BOOST_TEST(precice.getMeshDimensions(meshName2) == 2);
    BOOST_TEST(precice.getMeshDimensions(meshName3) == 2);

    std::vector<int> vertexIDs1;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName1, positions.at(i));
      vertexIDs1.push_back(vertexID);
    }
    std::vector<int> vertexIDs2;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName2, positions.at(i));
      vertexIDs2.push_back(vertexID);
    }
    std::vector<int> vertexIDs3;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName3, positions.at(i));
      vertexIDs3.push_back(vertexID);
    }

    precice.initialize();

    precice.writeData(meshName1, dataWriteID1, vertexIDs1, datas);
    precice.writeData(meshName2, dataWriteID2, vertexIDs2, datas);
    precice.writeData(meshName3, dataWriteID3, vertexIDs3, datas);

    if (precice.requiresWritingCheckpoint()) {
    }
    precice.advance(0.0001);
    if (precice.requiresReadingCheckpoint()) {
    }

    precice.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI

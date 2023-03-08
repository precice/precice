#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
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

  std::vector<Eigen::Vector2d> datas;
  Eigen::Vector2d              data;
  data << 1.0, 1.0;
  datas.push_back(data);
  data << 2.0, 2.0;
  datas.push_back(data);
  data << 3.0, 3.0;
  datas.push_back(data);
  data << 4.0, 5.0;
  datas.push_back(data);

  if (context.isNamed("SOLIDZ1") ||
      context.isNamed("SOLIDZ2") ||
      context.isNamed("SOLIDZ3")) {
    precice::SolverInterface precice(context.name, context.config(), 0, 1);
    BOOST_TEST(precice.getDimensions() == 2);

    std::string meshName, dataWriteID, dataReadID;
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

    std::vector<int> vertexIDs;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName, positions.at(i).data());
      vertexIDs.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i = 0; i < 4; i++) {
      precice.writeVectorData(meshName, dataWriteID, vertexIDs.at(i), datas.at(i).data());
    }

    if (precice.requiresWritingCheckpoint()) {
    }
    precice.advance(0.0001);
    if (precice.requiresReadingCheckpoint()) {
    }

    for (size_t i = 0; i < 4; i++) {
      precice.readVectorData(meshName, dataReadID, vertexIDs.at(i), datas.at(i).data());
    }

    BOOST_TEST(datas.at(0)(0) == 1.00000000000000002082e-03);
    BOOST_TEST(datas.at(0)(1) == 1.00000000000000002082e-03);
    BOOST_TEST(datas.at(1)(0) == 2.00000000000000000000e-03);
    BOOST_TEST(datas.at(1)(1) == 2.00000000000000002082e-03);
    BOOST_TEST(datas.at(2)(0) == 3.00000000000000006245e-03);
    BOOST_TEST(datas.at(2)(1) == 3.00000000000000006245e-03);
    BOOST_TEST(datas.at(3)(0) == 4.00000000000000008327e-03);
    BOOST_TEST(datas.at(3)(1) == 5.00000000000000010408e-03);

    precice.finalize();

  } else {
    BOOST_TEST(context.isNamed("NASTIN"));
    precice::SolverInterface precice("NASTIN", context.config(), 0, 1);
    BOOST_TEST(precice.getDimensions() == 2);
    auto meshName1    = "NASTIN_Mesh1";
    auto meshName2    = "NASTIN_Mesh2";
    auto meshName3    = "NASTIN_Mesh3";
    auto dataWriteID1 = "Forces1";
    auto dataWriteID2 = "Forces2";
    auto dataWriteID3 = "Forces3";

    std::vector<int> vertexIDs1;
    int              vertexID = -1;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName1, positions.at(i).data());
      vertexIDs1.push_back(vertexID);
    }
    std::vector<int> vertexIDs2;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName2, positions.at(i).data());
      vertexIDs2.push_back(vertexID);
    }
    std::vector<int> vertexIDs3;
    for (size_t i = 0; i < 4; i++) {
      vertexID = precice.setMeshVertex(meshName3, positions.at(i).data());
      vertexIDs3.push_back(vertexID);
    }

    precice.initialize();

    for (size_t i = 0; i < 4; i++) {
      precice.writeVectorData(meshName1, dataWriteID1, vertexIDs1.at(i), datas.at(i).data());
      precice.writeVectorData(meshName2, dataWriteID2, vertexIDs2.at(i), datas.at(i).data());
      precice.writeVectorData(meshName3, dataWriteID3, vertexIDs3.at(i), datas.at(i).data());
    }

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

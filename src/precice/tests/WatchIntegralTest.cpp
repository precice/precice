#include <Eigen/Core>
#include <algorithm>
#include <istream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include "../impl/WatchIntegral.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {
class Vertex;
} // namespace precice::mesh

using namespace precice;

namespace {
std::vector<double> readDoublesFromTXTFile(std::string const &filename, int skip = 0)
{
  std::ifstream is{filename};
  if (skip > 0) {
    std::string ignore;
    while (skip--) {
      is >> ignore;
    }
  }
  return {std::istream_iterator<double>{is}, std::istream_iterator<double>{}};
}
} // namespace

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(WatchIntegral)

BOOST_AUTO_TEST_CASE(ScalarDataNoConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh->createVertex(Eigen::Vector2d(1.0, 1.0));

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-noConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData NoConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 2);
    auto expected = std::vector<double>{
        0.0, 10.0,
        1.0, 14.0,
        2.0, 14.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataNoConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh->createVertex(Eigen::Vector2d(1.0, 1.0));

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4., 5., 6., 7., 8.}}));

  std::string fileName("precice-WatchIntegralTest-vectorData-noConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5., 6., 7., 8., 9.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData NoConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 16.0, 20.0,
        1.0, 20.0, 24.0,
        2.0, 20.0, 24.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataEdgeConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));

  mesh->createEdge(v1, v2);
  mesh->createEdge(v2, v3);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData EdgeConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 6.5, 3.0,
        1.0, 9.5, 3.0,
        2.0, 9.5, 3.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataEdgeConnectivityNoScale)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));

  mesh->createEdge(v1, v2);
  mesh->createEdge(v2, v3);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity-noScale.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, false);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData EdgeConnectivity NoScale")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 6.0, 3.0,
        1.0, 9.0, 3.0,
        2.0, 9.0, 3.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataEdgeConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));

  mesh->createEdge(v1, v2);
  mesh->createEdge(v2, v3);

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4., 5., 6.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5., 6., 7.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData EdgeConnectivity")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 4);
      auto expected = std::vector<double>{
          0.0, 10.0, 13.0, 3.0,
          1.0, 13.0, 56.0, 3.0,
          2.0, 13.0, 56.0, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataEdgeConnectivityNoScale)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));

  mesh->createEdge(v1, v2);
  mesh->createEdge(v2, v3);

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4., 5., 6.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity-noScale.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, false);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5., 6., 7.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData EdgeConnectivity NoScale")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 4);
      auto expected = std::vector<double>{
          0.0, 9.0, 12.0, 3.0,
          1.0, 12.0, 15.0, 3.0,
          2.0, 12.0, 15.0, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataFaceConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));

  mesh::Edge &e1 = mesh->createEdge(v1, v2);
  mesh::Edge &e2 = mesh->createEdge(v2, v3);
  mesh::Edge &e3 = mesh->createEdge(v3, v4);
  mesh::Edge &e4 = mesh->createEdge(v4, v1);
  mesh::Edge &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-faceConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData FaceConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 28.0, 12.0,
        1.0, 40.0, 12.0,
        2.0, 40.0, 12.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataFaceConnectivityNoScale)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));

  mesh::Edge &e1 = mesh->createEdge(v1, v2);
  mesh::Edge &e2 = mesh->createEdge(v2, v3);
  mesh::Edge &e3 = mesh->createEdge(v3, v4);
  mesh::Edge &e4 = mesh->createEdge(v4, v1);
  mesh::Edge &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4.}}));

  std::string fileName("precice-WatchIntegralTest-scalarData-faceConnectivity-noScale.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, false);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData FaceConnectivity NoScale")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 10.0, 12.0,
        1.0, 14.0, 12.0,
        2.0, 14.0, 12.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataFaceConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));

  mesh::Edge &e1 = mesh->createEdge(v1, v2);
  mesh::Edge &e2 = mesh->createEdge(v2, v3);
  mesh::Edge &e3 = mesh->createEdge(v3, v4);
  mesh::Edge &e4 = mesh->createEdge(v4, v1);
  mesh::Edge &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4., 5., 6., 7., 8.}}));

  std::string fileName("precice-WatchIntegralTest-vectorData-faceConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5., 6., 7., 8., 9.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData FaceConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 4);
    auto expected = std::vector<double>{
        0.0, 44.0, 56.0, 12.0,
        1.0, 56.0, 68.0, 12.0,
        2.0, 56.0, 68.0, 12.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataFaceConnectivityNoScale)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));

  mesh::Edge &e1 = mesh->createEdge(v1, v2);
  mesh::Edge &e2 = mesh->createEdge(v2, v3);
  mesh::Edge &e3 = mesh->createEdge(v3, v4);
  mesh::Edge &e4 = mesh->createEdge(v4, v1);
  mesh::Edge &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4., 5., 6., 7., 8.}}));

  std::string fileName("precice-WatchIntegralTest-vectorData-faceConnectivity-noScale.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, false);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4., 5., 6., 7., 8., 9.}}));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData FaceConnectivity NoScale")
  {
    auto result   = readDoublesFromTXTFile(fileName, 4);
    auto expected = std::vector<double>{
        0.0, 16.0, 20.0, 12.0,
        1.0, 20.0, 24.0, 12.0,
        2.0, 20.0, 24.0, 12.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(MeshChangeFaceConnectivity)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));

  mesh::Edge &e1 = mesh->createEdge(v1, v2);
  mesh::Edge &e2 = mesh->createEdge(v2, v3);
  mesh::Edge &e3 = mesh->createEdge(v3, v4);
  mesh::Edge &e4 = mesh->createEdge(v4, v1);
  mesh::Edge &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4.}}));

  std::string fileName("precice-WatchIntegralTest-meshChange-faceConnectivity.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    v2.setCoords(Eigen::Vector3d(3.0, -4.0, 0.0));

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral MeshChange FaceConnectivity")
  {
    auto result   = readDoublesFromTXTFile(fileName, 3);
    auto expected = std::vector<double>{
        0.0, 28.0, 12.0,
        1.0, 40.0, 18.0,
        2.0, 40.0, 18.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataNoConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));
  PtrData     doubleData = mesh->createData("DoubleData", 1, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  } else if (context.isRank(1)) {
    mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  } else if (context.isRank(2)) {
    mesh->createVertex(Eigen::Vector2d(2.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  } else if (context.isRank(3)) {
    mesh->createVertex(Eigen::Vector2d(3.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  }

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0}}));
  } else if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{3.0, 4.0}}));
  } else if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{5.0, 6.0}}));
  } else if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{7.0, 8.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-scalarData-noConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0}}));
    } else if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{4.0, 5.0}}));
    } else if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{6.0, 7.0}}));
    } else if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{8.0, 9.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData NoConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 2);
      auto expected = std::vector<double>{
          0.0, 36.0,
          1.0, 44.0,
          2.0, 44.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataNoConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));
  PtrData     doubleData = mesh->createData("DoubleData", 2, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  } else if (context.isRank(1)) {
    mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  } else if (context.isRank(2)) {
    mesh->createVertex(Eigen::Vector2d(2.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  } else if (context.isRank(3)) {
    mesh->createVertex(Eigen::Vector2d(3.0, 0.0));
    mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  }

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0, 3.0, 4.0}}));
  } else if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{5.0, 6.0, 7.0, 8.0}}));
  } else if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{9.0, 10.0, 11.0, 12.0}}));
  } else if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{13.0, 14.0, 15.0, 16.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-vectorData-noConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0, 4.0, 5.0}}));
    } else if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{6.0, 7.0, 8.0, 9.0}}));
    } else if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{10.0, 11.0, 12.0, 13.0}}));
    } else if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{14.0, 15.0, 16.0, 17.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData NoConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 3);
      auto expected = std::vector<double>{
          0.0, 64.0, 72.0,
          1.0, 72.0, 80.0,
          2.0, 72.0, 80.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataEdgeConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  if (utils::IntraComm::isPrimary()) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
    mesh->createEdge(v1, v2);
  }
  if (context.isRank(1)) {
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));
    mesh->createEdge(v3, v4);
  }
  if (context.isRank(2)) {
    mesh::Vertex &v5 = mesh->createVertex(Eigen::Vector2d(2.0, 1.0));
    mesh::Vertex &v6 = mesh->createVertex(Eigen::Vector2d(2.0, 2.0));
    mesh->createEdge(v5, v6);
  }
  if (context.isRank(3)) {
    mesh::Vertex &v7 = mesh->createVertex(Eigen::Vector2d(3.0, 1.0));
    mesh::Vertex &v8 = mesh->createVertex(Eigen::Vector2d(3.0, 3.0));
    mesh->createEdge(v7, v8);
  }

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0}}));
  }
  if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{3.0, 4.0}}));
  }
  if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{5.0, 6.0}}));
  }
  if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{7.0, 8.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0}}));
    }
    if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{4.0, 5.0}}));
    }
    if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{6.0, 7.0}}));
    }
    if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{8.0, 9.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData EdgeConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 3);
      auto expected = std::vector<double>{
          0.0, 25.5, 5.0,
          1.0, 30.5, 5.0,
          2.0, 30.5, 5.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataEdgeConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 2;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  if (utils::IntraComm::isPrimary()) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
    mesh->createEdge(v1, v2);
  }
  if (context.isRank(1)) {
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));
    mesh->createEdge(v3, v4);
  }
  if (context.isRank(2)) {
    mesh::Vertex &v5 = mesh->createVertex(Eigen::Vector2d(2.0, 1.0));
    mesh::Vertex &v6 = mesh->createVertex(Eigen::Vector2d(2.0, 2.0));
    mesh->createEdge(v5, v6);
  }
  if (context.isRank(3)) {
    mesh::Vertex &v7 = mesh->createVertex(Eigen::Vector2d(3.0, 1.0));
    mesh::Vertex &v8 = mesh->createVertex(Eigen::Vector2d(3.0, 3.0));
    mesh->createEdge(v7, v8);
  }

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0, 3.0, 4.0}}));
  }
  if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{5.0, 6.0, 7.0, 8.0}}));
  }
  if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{9.0, 10.0, 11.0, 12.0}}));
  }
  if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{13.0, 14.0, 15.0, 16.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-scalarData-edgeConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0, 4.0, 5.0}}));
    }
    if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{6.0, 7.0, 8.0, 9.0}}));
    }
    if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{10.0, 11.0, 12.0, 13.0}}));
    }
    if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{14.0, 15.0, 16.0, 17.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData0 DoubleData1  SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData EdgeConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 4);
      auto expected = std::vector<double>{
          0.0, 46.0, 51.0, 5.0,
          1.0, 51.0, 56.0, 5.0,
          2.0, 51.0, 56.0, 5.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScalarDataFaceConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  if (utils::IntraComm::isPrimary()) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
    mesh::Edge &  e1 = mesh->createEdge(v1, v2);
    mesh::Edge &  e2 = mesh->createEdge(v2, v3);
    mesh::Edge &  e5 = mesh->createEdge(v1, v3);
    mesh->createTriangle(e1, e2, e5);
  }
  if (context.isRank(1)) {
  }
  if (context.isRank(2)) {
  }
  if (context.isRank(3)) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));
    mesh::Edge &  e3 = mesh->createEdge(v3, v4);
    mesh::Edge &  e4 = mesh->createEdge(v4, v1);
    mesh::Edge &  e5 = mesh->createEdge(v1, v3);
    mesh->createTriangle(e3, e4, e5);
  }

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0, 3.0}}));
  }
  if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1));
  }
  if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1));
  }
  if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 3.0, 4.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-scalarData-faceConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0, 4.0}}));
    }
    if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1));
    }
    if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1));
    }
    if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 4.0, 5.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral ScalarData FaceConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 3);
      auto expected = std::vector<double>{
          0.0, 28.0, 12.0,
          1.0, 40.0, 12.0,
          2.0, 40.0, 12.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorDataFaceConnectivityParallel)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  using namespace mesh;
  // Setup geometry
  std::string name("rectangle");
  int         dimensions = 3;
  PtrMesh     mesh(new Mesh(name, dimensions, testing::nextMeshID()));

  if (utils::IntraComm::isPrimary()) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
    mesh::Edge &  e1 = mesh->createEdge(v1, v2);
    mesh::Edge &  e2 = mesh->createEdge(v2, v3);
    mesh::Edge &  e5 = mesh->createEdge(v1, v3);
    mesh->createTriangle(e1, e2, e5);
  }
  if (context.isRank(1)) {
  }
  if (context.isRank(2)) {
  }
  if (context.isRank(3)) {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 4.0, 0.0));
    mesh::Edge &  e3 = mesh->createEdge(v3, v4);
    mesh::Edge &  e4 = mesh->createEdge(v4, v1);
    mesh::Edge &  e5 = mesh->createEdge(v1, v3);
    mesh->createTriangle(e3, e4, e5);
  }

  PtrData doubleData = mesh->createData("DoubleData", 2, 0_dataID);

  if (utils::IntraComm::isPrimary()) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}}));
  }
  if (context.isRank(1)) {
    doubleData->setSampleAtTime(0, time::Sample(1));
  }
  if (context.isRank(2)) {
    doubleData->setSampleAtTime(0, time::Sample(1));
  }
  if (context.isRank(3)) {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0, 5.0, 6.0, 7.0, 8.0}}));
  }

  std::string fileName("precice-WatchIntegralTest-vectorData-faceConnectivity-parallel.log");

  {
    impl::WatchIntegral watchIntegral(mesh, fileName, true);
    watchIntegral.initialize();

    // Write output
    watchIntegral.exportIntegralData(0.0);

    // Change data (next timestep)
    if (utils::IntraComm::isPrimary()) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}));
    }
    if (context.isRank(1)) {
      doubleData->setSampleAtTime(1, time::Sample(1));
    }
    if (context.isRank(2)) {
      doubleData->setSampleAtTime(1, time::Sample(1));
    }
    if (context.isRank(3)) {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2.0, 3.0, 6.0, 7.0, 8.0, 9.0}}));
    }

    // Write output again
    watchIntegral.exportIntegralData(1.0);

    watchIntegral.exportIntegralData(2.0);
  }
  // File Format: Time  DoubleData SurfaceArea
  BOOST_TEST_CONTEXT("Validating WatchIntegral VectorData FaceConnectivity Parallel")
  {
    if (utils::IntraComm::isPrimary()) {
      auto result   = readDoublesFromTXTFile(fileName, 4);
      auto expected = std::vector<double>{
          0.0, 44.0, 56.0, 12.0,
          1.0, 56.0, 68.0, 12.0,
          2.0, 56.0, 68.0, 12.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END() // Precice

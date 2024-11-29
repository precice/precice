#include <Eigen/Core>
#include <algorithm>
#include <istream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include "../impl/WatchPoint.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {
class Vertex;
} // namespace precice::mesh

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PreciceTests)

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

void testWatchPoint(const TestContext & context,
                    bool                withEdge,
                    std::vector<double> watchPosition,
                    std::vector<double> expected)
{
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("rectangle");
  PtrMesh     mesh(new Mesh(name, 2, testing::nextMeshID()));

  if (context.size > 1) {
    if (context.isPrimary()) {
      mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
      mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
      if (withEdge) {
        mesh->createEdge(v1, v2);
      }
    } else {
      mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
      mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(0.0, 2.0));
      if (withEdge) {
        mesh->createEdge(v2, v3);
      }
    }
  } else {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(0.0, 2.0));
    if (withEdge) {
      mesh->createEdge(v1, v2);
      mesh->createEdge(v2, v3);
    }
  }

  using precice::testing::operator""_dataID;
  PtrData                 doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  PtrData                 vectorData = mesh->createData("VectorData", 2, 1_dataID);

  if (context.size > 1) {
    if (context.isPrimary()) {
      doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1.0, 2.0}}));
      vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{1.0, 2.0, 3.0, 4.0}}));
    } else {
      doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{2.0, 3.0}}));
      vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{3.0, 4.0, 5.0, 6.0}}));
    }
  } else {
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3.}}));
    vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{1., 2., 3., 4., 5., 6.}}));
  }

  std::string filename0("precice-WatchPointTest-timeseries-0.log");

  bool isWatchpointClosest = false;
  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector2d  pointToWatch0(watchPosition[0], watchPosition[1]);
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);

    // Initialize
    watchpoint0.initialize();
    isWatchpointClosest = watchpoint0.isClosest();
    // Write output
    watchpoint0.exportPointData(0.0);

    // Change data (next timestep)
    if (context.size > 1) {
      if (context.isPrimary()) {
        doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3.}}));
        vectorData->setSampleAtTime(1, time::Sample(2, Eigen::VectorXd{{2., 3., 4., 5.}}));
      } else {
        doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{3., 4.}}));
        vectorData->setSampleAtTime(1, time::Sample(2, Eigen::VectorXd{{4., 5., 6., 7.}}));
      }
    } else {
      doubleData->setSampleAtTime(1, time::Sample(1, Eigen::VectorXd{{2., 3., 4.}}));
      vectorData->setSampleAtTime(1, time::Sample(2, Eigen::VectorXd{{2., 3., 4., 5., 6., 7.}}));
    }

    // Write output again
    watchpoint0.exportPointData(1.0);

    // Write output again to check if the data stays the same
    watchpoint0.exportPointData(2.0);
  }

  // File Format: Time  Coordinate0  Coordinate1  DoubleData  VectorData0  VectorData1
  if (isWatchpointClosest) {
    BOOST_TEST_CONTEXT("Validating watchpoint0")
    {
      auto result = readDoublesFromTXTFile(filename0, 6);
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

BOOST_AUTO_TEST_SUITE(WatchPoint)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TimeSeriesNoEdgeSerialPoint1)
{
  PRECICE_TEST();
  bool withEdge           = false;
  auto watchPointPosition = std::vector<double>{-0.5, 0.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      1.0, 0.0, 1.0, 3.0, 4.0, 5.0,
      2.0, 0.0, 1.0, 3.0, 4.0, 5.0};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TimeSeriesWithEdgeSerialPoint1)
{
  PRECICE_TEST();
  bool withEdge           = true;
  auto watchPointPosition = std::vector<double>{-0.5, 0.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 0.6, 1.6, 2.2, 3.2,
      1.0, 0.0, 0.6, 2.6, 3.2, 4.2,
      2.0, 0.0, 0.6, 2.6, 3.2, 4.2};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TimeSeriesNoEdgeSerialPoint2)
{
  PRECICE_TEST();
  bool withEdge           = false;
  auto watchPointPosition = std::vector<double>{0.0, 1.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 2.0, 3.0, 5.0, 6.0,
      1.0, 0.0, 2.0, 4.0, 6.0, 7.0,
      2.0, 0.0, 2.0, 4.0, 6.0, 7.0};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TimeSeriesWithEdgeSerialPoint2)
{
  PRECICE_TEST();
  bool withEdge           = true;
  auto watchPointPosition = std::vector<double>{0.0, 1.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 1.6, 2.6, 4.2, 5.2,
      1.0, 0.0, 1.6, 3.6, 5.2, 6.2,
      2.0, 0.0, 1.6, 3.6, 5.2, 6.2};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TimeSeriesNoEdgeParallelPoint1)
{
  PRECICE_TEST();
  bool withEdge           = false;
  auto watchPointPosition = std::vector<double>{-0.5, 0.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      1.0, 0.0, 1.0, 3.0, 4.0, 5.0,
      2.0, 0.0, 1.0, 3.0, 4.0, 5.0};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TimeSeriesWithEdgeParallelPoint1)
{
  PRECICE_TEST();
  bool withEdge           = true;
  auto watchPointPosition = std::vector<double>{-0.5, 0.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 0.6, 1.6, 2.2, 3.2,
      1.0, 0.0, 0.6, 2.6, 3.2, 4.2,
      2.0, 0.0, 0.6, 2.6, 3.2, 4.2};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TimeSeriesNoEdgeParallelPoint2)
{
  PRECICE_TEST();
  bool withEdge           = false;
  auto watchPointPosition = std::vector<double>{0.0, 1.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 2.0, 3.0, 5.0, 6.0,
      1.0, 0.0, 2.0, 4.0, 6.0, 7.0,
      2.0, 0.0, 2.0, 4.0, 6.0, 7.0};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TimeSeriesWithEdgeParallelPoint2)
{
  PRECICE_TEST();
  bool withEdge           = true;
  auto watchPointPosition = std::vector<double>{0.0, 1.6};
  auto expected           = std::vector<double>{
      0.0, 0.0, 1.6, 2.6, 4.2, 5.2,
      1.0, 0.0, 1.6, 3.6, 5.2, 6.2,
      2.0, 0.0, 1.6, 3.6, 5.2, 6.2};
  testWatchPoint(context, withEdge, watchPointPosition, expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Reinitialize)
{
  PRECICE_TEST();
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("rectangle");
  PtrMesh     mesh(new Mesh(name, 2, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh->createEdge(v1, v2);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  PtrData vectorData = mesh->createData("VectorData", 2, 1_dataID);

  // v1, v2 carry data 1
  // v2 and later v3 carry data 2
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 1., 2.}}));
  vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{1., 1., 1., 1., 2., 2.}}));

  std::string filename0("precice-WatchPointTest-reinit-0.log");
  std::string filename1("precice-WatchPointTest-reinit-1.log");

  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector2d  pointToWatch0(0.1, 0.5);
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);
    Eigen::Vector2d  pointToWatch1(1.1, 0.5);
    impl::WatchPoint watchpoint1(pointToWatch1, mesh, filename1);

    // Initialize
    watchpoint0.initialize();
    watchpoint1.initialize();

    // Write output
    watchpoint0.exportPointData(0.0);
    watchpoint1.exportPointData(0.0);

    // Change Mesh and data
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
    mesh->createEdge(v3, v4);
    mesh->index().clear();
    doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 1., 2., 2.}}));
    vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{1., 1., 1., 1., 2., 2., 2., 2.}}));

    // Re-Initialize
    watchpoint0.initialize();
    watchpoint1.initialize();

    // Write output
    watchpoint0.exportPointData(1.0);
    watchpoint1.exportPointData(1.0);
  }

  // File Format: Time  Coordinate0  Coordinate1  DoubleData  VectorData0  VectorData1
  BOOST_TEST_CONTEXT("Validating watchpoint0")
  {
    auto result   = readDoublesFromTXTFile(filename0, 6);
    auto expected = std::vector<double>{
        0.0, 0.0, 0.5, 1.0, 1.0, 1.0,
        1.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }

  BOOST_TEST_CONTEXT("Validating watchpoint1")
  {
    auto result   = readDoublesFromTXTFile(filename1, 6);
    auto expected = std::vector<double>{
        0.0, 0.0, 0.5, 1.0, 1.0, 1.0,
        1.0, 1.0, 0.5, 2.0, 2.0, 2.0};
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VolumetricInterpolation2D)
{
  PRECICE_TEST();
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("triangle");
  PtrMesh     mesh(new Mesh(name, 2, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh->createEdge(v1, v2);
  mesh->createEdge(v2, v3);
  mesh->createEdge(v3, v1);
  mesh->createTriangle(v1, v2, v3);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);
  PtrData vectorData = mesh->createData("VectorData", 2, 1_dataID);

  // Data is (1,1,2) for the scalar, and same for each vector component.
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 1., 2., 2.}}));
  vectorData->setSampleAtTime(0, time::Sample(2, Eigen::VectorXd{{1., 1., 1., 1., 2., 2., 2., 2.}}));

  std::string filename0("precice-WatchPointTest-volumetric2d-0.log");
  std::string filename1("precice-WatchPointTest-volumetric2d-1.log");

  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector2d  pointToWatch0(0.1, 0.5);
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);
    Eigen::Vector2d  pointToWatch1(0.25, 0.25);
    impl::WatchPoint watchpoint1(pointToWatch1, mesh, filename1);

    // Initialize
    watchpoint0.initialize();
    watchpoint1.initialize();

    // Write output
    watchpoint0.exportPointData(0.0);
    watchpoint1.exportPointData(0.0);
  }

  // File Format: Time  Coordinate0  Coordinate1  DoubleData  VectorData0  VectorData1
  BOOST_TEST_CONTEXT("Validating watchpoint0")
  {
    auto result   = readDoublesFromTXTFile(filename0, 6);
    auto expected = std::vector<double>{
        0.0, 0.1, 0.5, 1.5, 1.5, 1.5};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }

  BOOST_TEST_CONTEXT("Validating watchpoint1")
  {
    auto result   = readDoublesFromTXTFile(filename1, 6);
    auto expected = std::vector<double>{
        0.0, 0.25, 0.25, 1.25, 1.25, 1.25};
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VolumetricInterpolation3D)
{
  PRECICE_TEST();
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("triangle");
  PtrMesh     mesh(new Mesh(name, 3, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  mesh->createTetrahedron(v1, v2, v3, v4);

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);

  // Data is (1,1,2) for the scalar, and same for each vector component.
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd{{1., 2., 3., 4.}}));

  std::string filename0("precice-WatchPointTest-volumetric3d-0.log");
  std::string filename1("precice-WatchPointTest-volumetric3d-1.log");

  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector3d  pointToWatch0(0.1, 0.5, 0.2);
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);
    Eigen::Vector3d  pointToWatch1(0.25, 0.25, 0.25);
    impl::WatchPoint watchpoint1(pointToWatch1, mesh, filename1);

    // Initialize
    watchpoint0.initialize();
    watchpoint1.initialize();

    // Write output
    watchpoint0.exportPointData(0.0);
    watchpoint1.exportPointData(0.0);
  }

  // File Format: Time  Coordinate0  Coordinate1 Coordinate2 DoubleData
  BOOST_TEST_CONTEXT("Validating watchpoint0")
  {
    auto result   = readDoublesFromTXTFile(filename0, 5);
    auto expected = std::vector<double>{
        0.0, 0.1, 0.5, 0.2, 2.7};
    BOOST_TEST(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
      BOOST_TEST_CONTEXT("entry index: " << i)
      {
        using testing::equals;
        BOOST_TEST(equals(result.at(i), expected.at(i)));
      }
    }
  }

  BOOST_TEST_CONTEXT("Validating watchpoint1")
  {
    auto result   = readDoublesFromTXTFile(filename1, 5);
    auto expected = std::vector<double>{
        0.0, 0.25, 0.25, 0.25, 2.5};
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

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(VolumetricParallel)
{
  PRECICE_TEST();
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("splitTetra");
  PtrMesh     mesh(new Mesh(name, 3, testing::nextMeshID()));

  switch (context.rank) {
  case 0: {
    mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  } break;
  case 1: {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
    mesh->createEdge(v1, v2);
  } break;
  case 2: {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
    mesh->createTriangle(v1, v2, v3);
  } break;

  case 3: {
    mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
    mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
    mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

    mesh->createTetrahedron(v1, v2, v3, v4);
  } break;
  }

  PtrData doubleData = mesh->createData("DoubleData", 1, 0_dataID);

  // Data is (1, 2, 3, 4) on the tetra, other ranks agree on their subset
  doubleData->setSampleAtTime(0, time::Sample(1, Eigen::VectorXd(context.rank + 1).setLinSpaced(1., context.rank + 1)));

  std::string filename0("precice-WatchPointTest-volumetricParallel-0.log");
  bool        isClosest;

  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector3d pointToWatch0(0.1, 0.2, 0.3);
    // Barycentric coordinates are (0.4, 0.1, 0.2, 0.3)
    // Thus expected value inside tetra is 0.4 + 0.1*2 + 0.2*3 + 0.3*4 = 2.4
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);

    watchpoint0.initialize();
    watchpoint0.exportPointData(0.0);
    isClosest = watchpoint0.isClosest();
  }

  // File Format: Time  Coordinate0  Coordinate1 Coordinate2 DoubleData
  // Closest rank does the check so we know the file is closed.
  if (isClosest) {
    BOOST_TEST_CONTEXT("Validating watchpoint0")
    {
      auto result   = readDoublesFromTXTFile(filename0, 5);
      auto expected = std::vector<double>{
          0.0, 0.1, 0.2, 0.3, 2.4};
      BOOST_TEST(result.size() = expected.size());
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

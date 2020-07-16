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

namespace precice {
namespace mesh {
class Vertex;
} // namespace mesh
} // namespace precice

using namespace precice;

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

BOOST_AUTO_TEST_SUITE(WatchPoint)

BOOST_AUTO_TEST_CASE(TimeSeries)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("rectangle");
  bool        flipNormals = false;
  PtrMesh     mesh(new Mesh(name, 2, flipNormals, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh->createEdge(v1, v2);

  PtrData doubleData   = mesh->createData("DoubleData", 1);
  PtrData vectorData   = mesh->createData("VectorData", 2);
  auto &  doubleValues = doubleData->values();
  auto &  vectorValues = vectorData->values();
  mesh->computeState();
  mesh->allocateDataValues();

  doubleValues(0) = 1.0;
  doubleValues(1) = 2.0;

  vectorValues(0) = 1.0;
  vectorValues(1) = 2.0;
  vectorValues(2) = 3.0;
  vectorValues(3) = 4.0;

  std::string filename0("precice-WatchPointTest-timeseries-0.log");
  std::string filename1("precice-WatchPointTest-timeseries-1.log");

  // this scope forces the filestreams to be closed
  {
    // Create watchpoints
    Eigen::Vector2d  pointToWatch0(1.0, 0.0); // always maps to v1
    impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);
    Eigen::Vector2d  pointToWatch1(-0.5, 0.5); // always maps to midpoint between v1 and v2
    impl::WatchPoint watchpoint1(pointToWatch1, mesh, filename1);

    // Initialize
    watchpoint0.initialize();
    watchpoint1.initialize();

    // Write output
    watchpoint0.exportPointData(0.0);
    watchpoint1.exportPointData(0.0);

    // Change data (next timestep)
    doubleValues(0) = 2.0;
    doubleValues(1) = 3.0;

    vectorValues(0) = 2.0;
    vectorValues(1) = 3.0;
    vectorValues(2) = 4.0;
    vectorValues(3) = 5.0;

    // Write output again
    watchpoint0.exportPointData(1.0);
    watchpoint1.exportPointData(1.0);

    // Write output again to check if the data stays the same
    watchpoint0.exportPointData(2.0);
    watchpoint1.exportPointData(2.0);
  }

  // File Format: Time  Coordinate0  Coordinate1  DoubleData  VectorData0  VectorData1
  BOOST_TEST_CONTEXT("Validating watchpoint0")
  {
    auto result   = readDoublesFromTXTFile(filename0, 6);
    auto expected = std::vector<double>{
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        2.0,
        1.0,
        0.0,
        0.0,
        2.0,
        2.0,
        3.0,
        2.0,
        0.0,
        0.0,
        2.0,
        2.0,
        3.0,
    };
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
        0.0, 0.0, 0.5, 1.5, 2.0, 3.0,
        1.0, 0.0, 0.5, 2.5, 3.0, 4.0,
        2.0, 0.0, 0.5, 2.5, 3.0, 4.0};
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

BOOST_AUTO_TEST_CASE(Reinitalize)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("rectangle");
  bool        flipNormals = false;
  PtrMesh     mesh(new Mesh(name, 2, flipNormals, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh->createEdge(v1, v2);

  PtrData doubleData   = mesh->createData("DoubleData", 1);
  PtrData vectorData   = mesh->createData("VectorData", 2);
  auto &  doubleValues = doubleData->values();
  auto &  vectorValues = vectorData->values();
  mesh->computeState();
  mesh->allocateDataValues();

  // v1, v2 carry data 1
  // v2 and later v3 carry data 2
  doubleValues.setConstant(1.0);
  doubleValues(2) = 2.0;
  vectorValues.setConstant(1.0);
  vectorValues(4) = 2.0;
  vectorValues(5) = 2.0;

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
    mesh->allocateDataValues();
    mesh->computeState();
    doubleValues.setConstant(1.0);
    doubleValues(2) = 2.0;
    doubleValues(3) = 2.0;
    vectorValues.setConstant(1.0);
    vectorValues(4) = 2.0;
    vectorValues(5) = 2.0;
    vectorValues(6) = 2.0;
    vectorValues(7) = 2.0;

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END() // Precice

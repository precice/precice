#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <string_view>
#include "logging/LogMacros.hpp"
#include "mapping/CoarseGrainingMapping.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/impl/MappingDataCache.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Sample.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(CoarseGrainingMapping)

namespace {

// Checks: outside radius -> 0, inside (<= radius) -> > 0
// Assumes 5x5x5 grid with unit spacing starting at (0,0,0)
// and linear indexing: idx = k*(Nx*Ny) + j*Nx + i
void check_support_5x5x5(const Eigen::MatrixXd &values,
                         double                 radius,
                         const Eigen::Vector3d &center)
{
  const int Nx = 5, Ny = 5, Nz = 5;
  const int N = Nx * Ny * Nz;

  BOOST_TEST(static_cast<int>(values.cols()) == N);

  const double r2 = radius * radius;

  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        const int idx = k * (Nx * Ny) + j * Nx + i;

        // grid point coordinate (unit spacing, origin at 0,0,0)
        Eigen::Vector3d p{static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)};
        const double    d = (p - center).norm();

        const Eigen::VectorXd v = values.col(idx);
        Eigen::VectorXd       zero(v.size());
        zero.setZero();

        if (d <= radius) {
          BOOST_TEST(v > zero, boost::test_tools::per_element());
        } else {
          BOOST_TEST(v == zero, boost::test_tools::per_element());
        }
      }
    }
  }
}
} // namespace

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Conservative3D)
{
  PRECICE_TEST();
  const int dimensions = 3;
  const int scalar     = 1;
  const int vector     = 3;

  // Create mesh to map from
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));

  const int    Nx = 5, Ny = 5, Nz = 5;
  const double spacing = 1.0;
  const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        Eigen::Vector3d p{x0 + spacing * i,
                          y0 + spacing * j,
                          z0 + spacing * k};
        outMesh->createVertex(p);
      }
    }
  }

  // Setup mapping with mapping coordinates and geometry used
  PtrMesh                                 inMesh = MeshConfiguration::getJustInTimeMappingMesh(dimensions);
  const double                            radius = 2.5;
  precice::mapping::CoarseGrainingMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  BOOST_TEST(mapping.isJustInTimeMapping() == true);

  mapping::impl::MappingDataCache dummyCache(1);

  // Constructs the index tree for the defined output mesh
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  Eigen::MatrixXd inCoords(dimensions, 1);
  inCoords.setConstant(0.0);
  Eigen::MatrixXd inData(scalar, 1);
  inData.setConstant(1.0);
  Eigen::MatrixXd outData(scalar, Nx * Ny * Nz);
  outData.setZero();
  mapping.mapConservativeAt(inCoords, inData, dummyCache, outData);

  // Test that we only have support within the function radius
  check_support_5x5x5(outData, radius, Eigen::Vector3d(0.0, 0.0, 0.0));
  const double zeroEvaluation = 0.133690152197192;
  BOOST_TEST(outData(0, 0) == zeroEvaluation);
  const double sqrtOneEvaluation = 0.06352956032410564;
  // evaluation with std::sqrt(1.0)
  BOOST_TEST(outData(0, 1) == sqrtOneEvaluation);

  // We mirror the contribution of the 0,0,0 point to the 1,0,0 point using the 2,0,0 point
  inCoords(0, 0) = 2.0;
  mapping.mapConservativeAt(inCoords, inData, dummyCache, outData);
  BOOST_TEST(outData(0, 1) == 2 * sqrtOneEvaluation);
  outData.setZero();

  // Testing in the middle of the domain
  inCoords.setConstant(2.5);
  mapping.mapConservativeAt(inCoords, inData, dummyCache, outData);
  BOOST_TEST(outData(0, 0) == 0.0);
  check_support_5x5x5(outData, radius, Eigen::Vector3d(2.5, 2.5, 2.5));

  // Testing vector data
  inData.resize(vector, 2);
  inData.row(0).setConstant(1);
  inData.row(1).setConstant(2);
  inData.row(2).setConstant(3);
  inCoords.resize(dimensions, 2);
  inCoords.setZero();
  outData.resize(vector, Nx * Ny * Nz);
  outData.setZero();

  mapping.mapConservativeAt(inCoords, inData, dummyCache, outData);
  check_support_5x5x5(outData, radius, Eigen::Vector3d(0.0, 0.0, 0.0));

  BOOST_TEST(outData(0, 0) == 2 * zeroEvaluation);
  BOOST_TEST(outData(1, 0) == 2 * 2 * zeroEvaluation);
  BOOST_TEST(outData(2, 0) == 2 * 3 * zeroEvaluation);

  BOOST_TEST(outData(1, 2) == 2 * outData(0, 2));
  BOOST_TEST(outData(2, 2) == 3 * outData(0, 2));
}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

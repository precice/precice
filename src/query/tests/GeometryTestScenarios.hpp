#ifndef PRECICE_QUERY_TESTS_GEOMETRYTESTSCENARIOS_HPP_
#define PRECICE_QUERY_TESTS_GEOMETRYTESTSCENARIOS_HPP_

#include "mesh/Mesh.hpp"
#include "utils/Dimensions.hpp"
#include <vector>

namespace precice {
namespace query {
namespace tests {

/**
 * @brief Provides test scenarios for geometry related functionality of solver interface.
 */
class GeometryTestScenarios
{
public:

  GeometryTestScenarios();

  /**
   * @brief Destructor. Destroys created scenarios.
   */
  ~GeometryTestScenarios();

  /**
   * @brief Test scenario for point queries.
   */
  struct PointQueryScenario {
    mesh::Mesh mesh;
    std::list<Eigen::VectorXd> queryCoords;
    std::list<double> validDistances;
    std::list<Eigen::VectorXd> validDistanceVectors;
    PointQueryScenario ( int dim );
  };

  /**
   * @brief (Creates) and returns point query test scenario.
   */
  const PointQueryScenario& pointQueryScenario ( int dim );

  /**
   * @brief Test scenario for position queries.
   */
  struct PositionQueryScenario {
    mesh::Mesh mesh;
    std::list<Eigen::VectorXd> queryCoords;
    std::list<int> validPositions;
    PositionQueryScenario ( int dim );
  };

  /**
   * @brief (Creates) and returns position query test scenario.
   */
  const PositionQueryScenario& positionQueryScenario ( int dim );

  /**
   * @brief Test scenario for voxel queries.
   */
  struct VoxelQueryScenario {
    mesh::Mesh mesh;
    std::list<Eigen::VectorXd> queryCenters;
    std::list<Eigen::VectorXd> queryHalflengths;
    std::list<bool> includeBoundaries;
    std::list<int> validPositions;
    VoxelQueryScenario ( int dim );
  };

  /**
   * @brief (Creates) and returns voxel query test scenario.
   */
  const VoxelQueryScenario& voxelQueryScenario ( int dim );

private:

  // @brief Test scenario for point queries.
  PointQueryScenario* _pointQueryScenario;

  // @brief Test scenario for position queries.
  PositionQueryScenario* _positionQueryScenario;

  // @brief Test scenarion for voxel queries.
  VoxelQueryScenario* _voxelQueryScenario;

  /**
   * @brief Creates a square (int 2D) or cube (in 3D) mesh.
   */
  void createMesh ( mesh::Mesh& mesh );

  /**
   * @brief Adds a query with results to the point query scenario.
   */
  void addToPointQueryScenario (
    const Eigen::VectorXd& queryCoord,
    double                 validDistance,
    const Eigen::VectorXd& validDistanceVector );
};

}}} // namespace precice, query, tests

#endif /* PRECICE_QUERY_TESTS_GEOMETRYTESTSCENARIOS_HPP_ */

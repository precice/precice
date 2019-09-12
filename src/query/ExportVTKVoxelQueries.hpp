#ifndef PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_
#define PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_

#include <Eigen/Core>
#include <string>
#include <vector>

namespace precice {
namespace mesh {
class Group;
class IVisitable;
} // namespace mesh
namespace query {

/**
 * @brief Visualizes results of FindVoxelContainedVisitor searches
 *
 * The results are visualized as voxels.
 *
 * @author Bernhard Gatzhammer
 */
class ExportVTKVoxelQueries {
public:
  void addQuery(
      const Eigen::VectorXd &voxelCenter,
      const Eigen::VectorXd &voxelHalflengths,
      int                    containedElementsCount);

  /**
    * @brief Writes performed queries' results into vtk file
    */
  void exportQueries(std::string filename);

  //   void exportVoxels ( std::string filename );

  /**
    * @brief Resets query results
    */
  void resetQueries();

private:
  // @brief Stores voxel centers queried
  std::vector<Eigen::VectorXd> _voxelCenters;

  // @brief Stores voxel halflengths queried
  std::vector<Eigen::VectorXd> _voxelHalflengths;

  // @brief Stores in queried voxels contained visitables
  std::vector<int> _containedElementCount;
};

} // namespace query
} // namespace precice

#endif /* PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_ */

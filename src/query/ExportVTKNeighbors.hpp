#pragma once

#include <string>
#include "FindClosest.hpp"

namespace precice {
namespace query {

/*
 * @brief Visualizes results of FindClosestOnMeshVisitor searches
 *
 * The visualization consists of lines connecting search points with
 * closest points on found visitables. The lines have assigned a value
 * distance, that has the value of the distance to the found visitables.
 *
 * @author Bernhard Gatzhammer
 */
class ExportVTKNeighbors {
public:
  /**
    * @brief Adds a neighborhood relation already computed.
    *
    * @param[in] queryPoint Point whose next neighbor was determined.
    * @param[in] closestNeighbor ClosestNeighbor of queryPoint.
    */
  void addNeighbors(
      const Eigen::VectorXd &queryPoint,
      const ClosestElement & closestNeighbor);

  /**
    * @brief Visualizes results of performed searches into vtk file
    */
  void exportNeighbors(const std::string &filename);

  /**
    * @brief Resets results of performed searches
    */
  void resetElements()
  {
    _neighbors.clear();
  };

private:
  /// Results of performed searches
  std::vector<std::pair<Eigen::VectorXd, query::ClosestElement>> _neighbors;
};

} // namespace query
} // namespace precice

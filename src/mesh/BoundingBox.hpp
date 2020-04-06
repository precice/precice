#pragma once
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace mesh {

/// Bounding box class for mesh and partions
class BoundingBox {

public:
  /***
  * 
  * @brief Constructor.
  * 
  * @param[in] bounds Vector of minimum and maximum points in each dimension
  * @param[in] safetyFactor Factor which enlarges the bounding box
  * 
  */
  BoundingBox(std::vector<double> bounds, double safetyFactor);
  ~BoundingBox();

  /// Bounding box factory function
  static BoundingBox createFromData(std::vector<double> bounds, double safetyFactor);

  /***
  * 
  * @brief Set bounds manually
  * 
  * @param[in] dimension Dimension in which the bounds are given
  * @param[in] min Minimum bound
  * @param[in] max Maximum bound
  * 
  */
  void setBounds(int dimension, double min, double max);

  /// Merges the bounding box with given bounding box, also enlage by _safetyFactor
  void mergeBoundingBoxes(const BoundingBox &otherBB, double safety);

  /// Checks if vertex in contained in _bb
  bool isVertexInBB(const mesh::Vertex &vertex);

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB);

  friend std::ostream &operator<<(std::ostream &out, const BoundingBox &bb);

private:
  logging::Logger _log{"mesh::BoundingBox"};

  BoundingBox(const BoundingBox &bb);

  int    _dimensions;
  double _safetyFactor;

  /// Whether this bounding box is prepared or not
  bool _prepared;

  /// Container of min and max points in each dimension
  std::array<double, 6> _bounds;
};

} // namespace mesh
} // namespace precice

#pragma once
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"
#include <vector>

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
  BoundingBox(std::vector<double> bounds);
  BoundingBox(int dimension);
  BoundingBox();
  ~BoundingBox();
  BoundingBox(const BoundingBox &bb);
  bool operator==(const BoundingBox& otherBB) const;

  /// Bounding box factory function
  static BoundingBox createFromData(std::vector<double> bounds);

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
  void setMin(int dimension, double min);
  void setMax(int dimension, double max);

  void setSafetyFactor(double safetyFactor);

  /// Merges the bounding box with given bounding box, also enlage by _safetyFactor
  bool mergeBoundingBoxes(const BoundingBox &otherBB);

  /// Expand bounding box
  void expandTo(const Vertex& vertices);

  /// Checks if vertex in contained in _bb
  bool isVertexInBB(const Vertex &vertex);

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB);

  /// Return the value for given dimension and bound type
  double getData(int dimension, int type) const;

  const double* data() const;
  std::vector<double> dataVector();

  bool empty();

  int getDimension() const;

  friend std::ostream &operator<<(std::ostream &out, const BoundingBox &bb);

private:
  logging::Logger _log{"mesh::BoundingBox"};

  int    _dimensions;
  double _safetyFactor;

  /// Whether this bounding box is prepared or not
  bool _prepared{false};

  /// Container of min and max points in each dimension
  std::vector<double> _bounds;
};

} // namespace mesh
} // namespace precice

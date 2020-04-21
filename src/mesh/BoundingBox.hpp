#pragma once
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"
#include <vector>

namespace precice {
namespace mesh {

/// An axis-aligned bounding box around a (partition of a) mesh
class BoundingBox {

public:
  /***
  * 
  * @brief Constructor.
  * 
  * @param[in] bounds Vector of minimum and maximum points in each dimension
  * 
  */
  BoundingBox(std::vector<double> bounds);
  
  /***
  * 
  * @brief Constructor.
  * 
  * @param[in] dimension Dimension of the bounding box
  * 
  */
  BoundingBox(int dimension);

  /// Copy Constructor
  BoundingBox(const BoundingBox &bb);

  /// Comparison Operator
  bool operator==(const BoundingBox& otherBB) const;

  /// Bounding box factory function
  static BoundingBox createFromData(std::vector<double> bounds);

  /// Expand bounding box using another bounding box
  void expandTo(const BoundingBox &otherBB);

  /// Expand bounding box using vertices
  void expandTo(const Vertex& vertices);

  /// Expand bounding box using value
  void expandTo(double value);

  /// Increase the size of bounding box by safety margin
  void addSafetyMargin(double safetyFactor);

  /// Checks if vertex in contained in _bb
  bool isVertexInBB(const Vertex &vertex) const;

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB);

  /**
   * @brief Returns the Center Of Gravity of the mesh
   *
   * Returns a vector of doubles, size d, each dimension computed as
   * cog =  (max - min) / 2 + min
   */
  std::vector<double> getCOG() const;

  /// Calculate the area of bounding box
  double getArea(std::vector<bool> deadAxis);

  /// Return data as std::vector
  const std::vector<double> dataVector() const;

  /// Getter dimension of the bounding box
  int getDimension() const;

  /// Output operator for easy logging
  friend std::ostream &operator<<(std::ostream &out, const BoundingBox &bb);

private:
  logging::Logger _log{"mesh::BoundingBox"};

  /// Number of dimensions
  int    _dimensions;

  /// Container of min and max points in each dimension
  std::vector<double> _bounds;
};

} // namespace mesh
} // namespace precice

#pragma once
#include <Eigen/Core>
#include <iosfwd>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace logging {
class Logger;
} // namespace logging

namespace mesh {
class Vertex;

/// An axis-aligned bounding box around a (partition of a) mesh
class BoundingBox {

public:
  /***
  *
  * @brief Constructor.
  *
  * @param[in] dimension Dimension of the bounding box
  *
  */
  explicit BoundingBox(int dimension);

  /***
  *
  * @brief Constructor.
  *
  * @param[in] bounds Min-max values of the bounding box in each dimesion
  *
  */
  explicit BoundingBox(std::vector<double> bounds);

  /// Special Members
  BoundingBox(const BoundingBox &) = default;
  BoundingBox(BoundingBox &&)      = default;
  BoundingBox &operator=(const BoundingBox &bb) = default;
  BoundingBox &operator=(BoundingBox &&bb) = default;

  /// Comparison Operator
  bool operator==(const BoundingBox &otherBB) const;

  /// Check whether the bounding box is at default state or not
  bool empty() const;

  /// Expand bounding box using another bounding box
  void expandBy(const BoundingBox &otherBB);

  /// Expand bounding box using vertices
  void expandBy(const Vertex &vertices);

  /// Expand bounding box using value
  void expandBy(double value);

  /// Increase the size of bounding box by safety margin
  void scaleBy(double safetyFactor);

  /// Checks if vertex in contained in _bb
  bool contains(const Vertex &vertex) const;

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB);

  /**
   * @brief Returns the Center Of Gravity of the mesh
   *
   * Returns a vector of doubles, size d, each dimension computed as
   * cog =  (max - min) / 2 + min
   */
  Eigen::VectorXd center() const;

  /// the min corner of the bounding box
  Eigen::VectorXd minCorner() const;

  /// the max corner of the bounding box
  Eigen::VectorXd maxCorner() const;

  /// Calculate the area of bounding box
  double getArea(std::vector<bool> deadAxis);

  /// Return data as std::vector
  const std::vector<double> dataVector() const;

  /// Getter dimension of the bounding box
  int getDimension() const;

  /// Print bounds of bounding box, output operator overload
  void print(std::ostream &out) const;

private:
  static logging::Logger _log;

  /// Number of dimensions
  int _dimensions;

  /// Container of min and max points in each dimension
  std::vector<double> _bounds;
};

std::ostream &operator<<(std::ostream &, const BoundingBox &);

} // namespace mesh
} // namespace precice

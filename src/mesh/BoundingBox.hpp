#pragma once
#include <Eigen/Core>
#include <fmt/ostream.h>
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
   * @brief Constructor
   *
   * @param[in] dimension Dimension of the bounding box
   *
   */
  explicit BoundingBox(int dimension);

  /***
   *
   * @brief Constructor.
   *
   * @param[in] bounds Min-max values of the bounding box in each dimension (x_min, x_max, y_min, y_max, z_min, z_max)
   *
   */
  explicit BoundingBox(std::vector<double> bounds);

  /**
   * @brief Constructor accepting Eigen::VectorXd
   *
   * @param boundMin: vertex at minCorner
   * @param boundMax: vertex at maxCorner
   */
  explicit BoundingBox(Eigen::VectorXd boundMin, Eigen::VectorXd boundMax);

  /// Special Members
  BoundingBox(const BoundingBox &)              = default;
  BoundingBox(BoundingBox &&)                   = default;
  BoundingBox &operator=(const BoundingBox &bb) = default;
  BoundingBox &operator=(BoundingBox &&bb)      = default;

  /// Comparison Operator
  bool operator==(const BoundingBox &otherBB) const;

  /// Check if every dimension's length is equal to zero
  bool empty() const;

  /// Check whether the bounding box is at default state or not
  /// all the values of _boundMin = std::numeric_limits<double>::max()
  /// all the values of _boundMax = std::numeric_limits<double>::lowest()
  bool isDefault() const;

  /// Expand bounding box using another bounding box
  void expandBy(const BoundingBox &otherBB);

  /// Expand bounding box using vertices
  void expandBy(const Vertex &vertices);

  /// Expand bounding box using a double value in all dimensions
  /// Using this method, BoundingBox should not be in the default state; otherwise, it will make ReceivedPartitionTest fail because of incorrect connectionMapSize
  void expandBy(double value);

  /// Increase the size of bounding box by safety margin
  void scaleBy(double safetyFactor);

  /// Checks if vertex in contained in _bb
  bool contains(const Vertex &vertex) const;

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB) const;

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

  /// returns the edge length of a specific axis
  double getEdgeLength(int axis) const;

  /// returns the maximum length of the bounding box in any dimension
  double longestEdgeLength() const;

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

  /// Number of dimensions (2 or 3)
  int _dimensions;

  /// _boundMin defining the bounding box, with minimum coordinates in each direction
  /// (x_min, y_min) when _dimensions=2;
  /// (x_min, y_min, z_min) when _dimensions=3;
  Eigen::VectorXd _boundMin;

  /// _boundMax defining the bounding box, with maximum coordinates in each direction
  /// (x_max, y_max) when _dimensions=2;
  /// (x_max, y_max, z_max) when _dimensions=3;
  Eigen::VectorXd _boundMax;
};

std::ostream &operator<<(std::ostream &, const BoundingBox &);

} // namespace mesh
} // namespace precice

template <>
struct fmt::formatter<precice::mesh::BoundingBox> : ostream_formatter {
};

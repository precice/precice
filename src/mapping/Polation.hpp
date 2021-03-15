#pragma once

#include <iosfwd>
#include <vector>
#include "Eigen/Core"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace mapping {

/// Struct that contains weight and index of a vertex
struct WeightedElement {
  int    vertexID;
  double weight;
};

/**
 * @brief Calculates the barycentric coordinates of a coordinate on the given vertex/edge/triangle and stores the corresponding weights
 * If all barycentric coordinates are positive, the operation is interpolation. If not, it is an extrapolation.
 */
class Polation {
public:
  /// Calculate projection to a vertex. Weight is always 1.0
  Polation(const mesh::Vertex &element);

  /// Calculate projection to an edge
  Polation(const Eigen::VectorXd &location, const mesh::Edge &element);

  /// Calculate projection to a triangle
  Polation(const Eigen::VectorXd &location, const mesh::Triangle &element);

  /// Get the weights and indices of the calculated interpolation
  const std::vector<WeightedElement> &getWeightedElements() const;

  /// Check whether all the weights are positive, which means it is interpolation
  bool isInterpolation() const;

private:
  std::vector<WeightedElement> _weightedElements;
};

/// Make the WeightedElement printable
std::ostream &operator<<(std::ostream &os, const WeightedElement &w);

/// Make the Polation class printable
std::ostream &operator<<(std::ostream &os, const Polation &p);

} // namespace mapping
} // namespace precice
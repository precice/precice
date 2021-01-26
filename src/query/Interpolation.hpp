#pragma once

#include "math/barycenter.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace query {

/**
 * @brief Weighting and reference to target element for a value to interpolate
 */
struct InterpolationElement {
  const mesh::Vertex *element = nullptr;
  double              weight  = 0.0;

  InterpolationElement() = default;
  InterpolationElement(const mesh::Vertex *element_, double weight_)
      : element(element_), weight(weight_) {}
  InterpolationElement(const mesh::Vertex &element_, double weight_)
      : element(&element_), weight(weight_) {}
};

/// Generates the InterpolationElements for directly projecting a Vertex on another Vertex
std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex &location,
    const mesh::Vertex &element);

/// Generates the InterpolationElements for projecting a Vertex on an Edge
std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex &location,
    const mesh::Edge &  element);

/// Generates the InterpolationElements for projecting a Vertex on a Triangle
std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex &  location,
    const mesh::Triangle &element);

} // namespace query
} // namespace precice
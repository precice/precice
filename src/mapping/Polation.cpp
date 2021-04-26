#include "mapping/Polation.hpp"
#include "math/barycenter.hpp"

namespace precice {
namespace mapping {

Polation::Polation(const mesh::Vertex &element)
{
  _weightedElements.emplace_back(WeightedElement{element.getID(), 1.0});
}

Polation::Polation(const Eigen::VectorXd &location, const mesh::Edge &element)
{
  PRECICE_ASSERT(location.size() == element.getDimensions(), location.size(), element.getDimensions());
  const auto &A = element.vertex(0);
  const auto &B = element.vertex(1);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForEdge(
                           A.getCoords(),
                           B.getCoords(),
                           element.computeNormal(),
                           location)
                           .barycentricCoords;

  _weightedElements.emplace_back(WeightedElement{A.getID(), bcoords(0)});
  _weightedElements.emplace_back(WeightedElement{B.getID(), bcoords(1)});
}

Polation::Polation(const Eigen::VectorXd &location, const mesh::Triangle &element)
{
  PRECICE_ASSERT(location.size() == element.getDimensions(), location.size(), element.getDimensions());
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);
  auto &C = element.vertex(2);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForTriangle(
                           A.getCoords(),
                           B.getCoords(),
                           C.getCoords(),
                           element.computeNormal(),
                           location)
                           .barycentricCoords;

  _weightedElements.emplace_back(WeightedElement{A.getID(), bcoords(0)});
  _weightedElements.emplace_back(WeightedElement{B.getID(), bcoords(1)});
  _weightedElements.emplace_back(WeightedElement{C.getID(), bcoords(2)});
}

const std::vector<WeightedElement> &Polation::getWeightedElements() const
{
  return _weightedElements;
}

bool Polation::isInterpolation() const
{
  return std::all_of(_weightedElements.begin(), _weightedElements.end(), [](const mapping::WeightedElement &elem) { return elem.weight >= 0.0; });
}

std::ostream &operator<<(std::ostream &os, const WeightedElement &w)
{
  return os << "(Vertex ID: " << w.vertexID << ", Weight: " << w.weight << ")";
}

std::ostream &operator<<(std::ostream &os, const Polation &p)
{
  os << "Polation: ";
  for (const auto &elem : p.getWeightedElements()) {
    os << elem;
  }
  return os;
}

} // namespace mapping
} // namespace precice

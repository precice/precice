#include <Eigen/Core>

#include "mapping/Polation.hpp"
#include "math/barycenter.hpp"
#include "math/differences.hpp"

namespace precice::mapping {

Polation::Polation(const Eigen::VectorXd &location, const mesh::Vertex &element)
{
  _weightedElements = {WeightedElement{element.getID(), 1.0}};
  // The projection in this case is simply the nearest point.
  _distance = (location - element.getCoords()).norm();
}

Polation::Polation(const Eigen::VectorXd &location, const mesh::Edge &element)
{
  PRECICE_ASSERT(location.size() == element.getDimensions(), location.size(), element.getDimensions());
  const auto &A = element.vertex(0);
  const auto &B = element.vertex(1);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForEdge(
      A.getCoords(),
      B.getCoords(),
      location);

  _weightedElements = {WeightedElement{A.getID(), bcoords(0)},
                       WeightedElement{B.getID(), bcoords(1)}};

  Eigen::VectorXd projection = A.getCoords() * bcoords(0) +
                               B.getCoords() * bcoords(1);
  _distance = (location - projection).norm();
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
      location);

  _weightedElements = {WeightedElement{A.getID(), bcoords(0)},
                       WeightedElement{B.getID(), bcoords(1)},
                       WeightedElement{C.getID(), bcoords(2)}};

  Eigen::VectorXd projection = A.getCoords() * bcoords(0) +
                               B.getCoords() * bcoords(1) +
                               C.getCoords() * bcoords(2);
  _distance = (location - projection).norm();
}

Polation::Polation(const Eigen::VectorXd &location, const mesh::Tetrahedron &element)
{
  PRECICE_ASSERT(location.size() == element.getDimensions(), location.size(), element.getDimensions());
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);
  auto &C = element.vertex(2);
  auto &D = element.vertex(3);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForTetrahedron(
      A.getCoords(),
      B.getCoords(),
      C.getCoords(),
      D.getCoords(),
      location);

  _weightedElements = {WeightedElement{A.getID(), bcoords(0)},
                       WeightedElement{B.getID(), bcoords(1)},
                       WeightedElement{C.getID(), bcoords(2)},
                       WeightedElement{D.getID(), bcoords(3)}};

  // There is no projection happening, so the distance is always 0.
  _distance = 0.0;
}

const boost::container::static_vector<WeightedElement, 4> &Polation::getWeightedElements() const
{
  return _weightedElements;
}

std::size_t Polation::nElements() const
{
  return _weightedElements.size();
}

bool Polation::isInterpolation() const
{
  return std::all_of(_weightedElements.begin(), _weightedElements.end(), [](const mapping::WeightedElement &elem) { return precice::math::greaterEquals(elem.weight, 0.0); });
}

double Polation::distance() const
{
  return _distance;
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

} // namespace precice::mapping

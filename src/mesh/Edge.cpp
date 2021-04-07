#include "Edge.hpp"
#include <Eigen/Core>
#include <algorithm>
#include "math/differences.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

Edge::Edge(
    Vertex &vertexOne,
    Vertex &vertexTwo,
    int     id)
    : _vertices({&vertexOne, &vertexTwo}),
      _id(id),
      _normal(Eigen::VectorXd::Constant(vertexOne.getDimensions(), 0.0))
{
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexTwo.getDimensions(),
                 vertexOne.getDimensions(), vertexTwo.getDimensions());
}

int Edge::getID() const
{
  return _id;
}

double Edge::getLength() const
{
  double length = (_vertices[1]->getCoords() - _vertices[0]->getCoords()).norm();
  return length;
}

const Eigen::VectorXd Edge::computeNormal(bool flip)
{
  // Compute normal
  Eigen::VectorXd edgeVector = vertex(1).getCoords() - vertex(0).getCoords();
  Eigen::VectorXd normal     = Eigen::Vector2d(-edgeVector[1], edgeVector[0]);
  if (not flip) {
    normal *= -1.0; // Invert direction if counterclockwise
  }
  _normal = normal.normalized(); // Scale normal vector to length 1

  return normal * getEnclosingRadius() * 2.0; // Weight by length
}

const Eigen::VectorXd Edge::getCenter() const
{
  return 0.5 * (_vertices[0]->getCoords() + _vertices[1]->getCoords());
}

double Edge::getEnclosingRadius() const
{
  return (_vertices[0]->getCoords() - getCenter()).norm();
}

bool Edge::connectedTo(const Edge &other) const
{
  return _vertices[0] == other._vertices[0] || _vertices[0] == other._vertices[1] || _vertices[1] == other._vertices[0] || _vertices[1] == other._vertices[1];
}

bool Edge::operator==(const Edge &other) const
{
  return math::equals(_normal, other._normal) &&
         std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin(),
                             [](const Vertex *a, const Vertex *b) { return *a == *b; });
}
bool Edge::operator!=(const Edge &other) const
{
  return !(*this == other);
}

std::ostream &operator<<(std::ostream &stream, const Edge &edge)
{
  using utils::eigenio::wkt;
  return stream << "LINESTRING ("
                << edge.vertex(0).getCoords().transpose().format(wkt())
                << ", "
                << edge.vertex(1).getCoords().transpose().format(wkt())
                << ')';
}

} // namespace mesh
} // namespace precice

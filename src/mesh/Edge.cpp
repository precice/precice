#include <Eigen/Core>
#include <algorithm>

#include "Edge.hpp"
#include "math/differences.hpp"
#include "precice/types.hpp"
#include "utils/EigenIO.hpp"

namespace precice::mesh {

Edge::Edge(
    Vertex &vertexOne,
    Vertex &vertexTwo)
    : _vertices({&vertexOne, &vertexTwo})
{
  if (*_vertices[1] < *_vertices[0]) {
    std::swap(_vertices[0], _vertices[1]);
  }
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexTwo.getDimensions(),
                 vertexOne.getDimensions(), vertexTwo.getDimensions());
}

double Edge::getLength() const
{
  double length = (_vertices[1]->getCoords() - _vertices[0]->getCoords()).norm();
  return length;
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
  return std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin(),
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

} // namespace precice::mesh

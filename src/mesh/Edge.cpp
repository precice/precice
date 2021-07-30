#include <Eigen/Core>
#include <algorithm>

#include "Edge.hpp"
#include "math/differences.hpp"
#include "precice/types.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

Edge::Edge(
    Vertex &vertexOne,
    Vertex &vertexTwo,
    EdgeID  id)
    : _vertices({&vertexOne, &vertexTwo}),
      _id(id)
{
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexTwo.getDimensions(),
                 vertexOne.getDimensions(), vertexTwo.getDimensions());
}

EdgeID Edge::getID() const
{
  return _id;
}

double Edge::getLength() const
{
  double length = (_vertices[1]->getCoords() - _vertices[0]->getCoords()).norm();
  return length;
}

Eigen::VectorXd Edge::computeNormal() const
{
  // Compute normal
  Eigen::VectorXd edgeVector = vertex(1).getCoords() - vertex(0).getCoords();

  // In 3D, use a normal on the plane z=0
  Eigen::VectorXd normal = Eigen::VectorXd::Zero(edgeVector.size());
  normal[0]              = -edgeVector[1];
  normal[1]              = edgeVector[0];

  // Handle vectors along the z axis
  if (normal.size() == 3 && normal.isApproxToConstant(0.0)) {
    normal[0] = 1.0;
  }

  return normal.normalized();
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

} // namespace mesh
} // namespace precice

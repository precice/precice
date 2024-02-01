#include "Tetrahedron.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <boost/concept/assert.hpp>
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenIO.hpp"

namespace precice::mesh {

Tetrahedron::Tetrahedron(
    Vertex &vertexOne,
    Vertex &vertexTwo,
    Vertex &vertexThree,
    Vertex &vertexFour)
    : _vertices({&vertexOne, &vertexTwo, &vertexThree, &vertexFour})
{
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexTwo.getDimensions(),
                 vertexOne.getDimensions(), vertexTwo.getDimensions());
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexThree.getDimensions(),
                 vertexOne.getDimensions(), vertexThree.getDimensions());
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexFour.getDimensions(),
                 vertexOne.getDimensions(), vertexFour.getDimensions());
  PRECICE_ASSERT(getDimensions() == 3, getDimensions());

  PRECICE_ASSERT(
      (&vertexOne != &vertexTwo) &&
          (&vertexOne != &vertexThree) &&
          (&vertexOne != &vertexFour) &&
          (&vertexTwo != &vertexThree) &&
          (&vertexTwo != &vertexFour) &&
          (&vertexThree != &vertexFour),
      "Tetrahedron vertices are not unique!");

  std::sort(_vertices.begin(), _vertices.end(),
            [](const Vertex *lhs, const Vertex *rhs) { return *lhs < *rhs; });
}

double Tetrahedron::getVolume() const
{
  return math::geometry::tetraVolume(vertex(0).getCoords(), vertex(1).getCoords(), vertex(2).getCoords(), vertex(3).getCoords());
}

int Tetrahedron::getDimensions() const
{
  return _vertices[0]->getDimensions();
}

const Eigen::VectorXd Tetrahedron::getCenter() const
{
  return (vertex(0).getCoords() + vertex(1).getCoords() + vertex(2).getCoords() + vertex(3).getCoords()) / 4.0;
}

double Tetrahedron::getEnclosingRadius() const
{
  auto center = getCenter();
  return std::max({(center - vertex(0).getCoords()).norm(),
                   (center - vertex(1).getCoords()).norm(),
                   (center - vertex(2).getCoords()).norm(),
                   (center - vertex(3).getCoords()).norm()});
}

bool Tetrahedron::operator==(const Tetrahedron &other) const
{
  return std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin(),
                             [](const Vertex *v1, const Vertex *v2) { return *v1 == *v2; });
}

bool Tetrahedron::operator!=(const Tetrahedron &other) const
{
  return !(*this == other);
}

std::ostream &operator<<(std::ostream &os, const Tetrahedron &t)
{
  // Show 6 edges: 0-1, 0-2, 0-3, 1-2, 1-3, 2_3
  using utils::eigenio::wkt;
  const auto &v0 = t.vertex(0);
  const auto &v1 = t.vertex(1);
  const auto &v2 = t.vertex(2);
  const auto &v3 = t.vertex(3);

  return os << "MULTILINESTRING ("
            << "(" << v0.getCoords().transpose().format(wkt()) << ", " << v1.getCoords().transpose().format(wkt()) << "), "
            << "(" << v0.getCoords().transpose().format(wkt()) << ", " << v2.getCoords().transpose().format(wkt()) << "), "
            << "(" << v0.getCoords().transpose().format(wkt()) << ", " << v3.getCoords().transpose().format(wkt()) << "), "
            << "(" << v1.getCoords().transpose().format(wkt()) << ", " << v2.getCoords().transpose().format(wkt()) << "), "
            << "(" << v1.getCoords().transpose().format(wkt()) << ", " << v3.getCoords().transpose().format(wkt()) << "), "
            << "(" << v2.getCoords().transpose().format(wkt()) << ", " << v3.getCoords().transpose().format(wkt()) << ")"
            << ")";
}

} // namespace precice::mesh

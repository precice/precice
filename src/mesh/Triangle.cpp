#include "Triangle.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <boost/concept/assert.hpp>
#include <boost/range/concepts.hpp>
#include <iterator>
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenIO.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Triangle::iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Triangle::const_iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<Triangle>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<const Triangle>) );

Triangle::Triangle(
    Edge &edgeOne,
    Edge &edgeTwo,
    Edge &edgeThree)
{
  PRECICE_ASSERT(edgeOne.getDimensions() == edgeTwo.getDimensions(),
                 edgeOne.getDimensions(), edgeTwo.getDimensions());
  PRECICE_ASSERT(edgeTwo.getDimensions() == edgeThree.getDimensions(),
                 edgeTwo.getDimensions(), edgeThree.getDimensions());

  PRECICE_ASSERT(edgeOne.connectedTo(edgeTwo), "Edge one and two are not connected.");
  PRECICE_ASSERT(edgeOne.connectedTo(edgeThree), "Edge one and three are not connected.");
  PRECICE_ASSERT(edgeTwo.connectedTo(edgeThree), "Edge two and three are not connected.");

  // Pick first 2 vertices from first edge
  Vertex &v0 = edgeOne.vertex(0);
  Vertex &v1 = edgeOne.vertex(1);

  // Determine the third vertex using the second edge
  Vertex *v2 = nullptr;
  if ((v0 == edgeTwo.vertex(0)) || (v1 == edgeTwo.vertex(0))) {
    v2 = &edgeTwo.vertex(1);
  } else if ((v0 == edgeTwo.vertex(1)) || (v1 == edgeTwo.vertex(1))) {
    v2 = &edgeTwo.vertex(0);
  } else {
    PRECICE_UNREACHABLE("Edges don't form a triangle");
  }

  _vertices = {&v0, &v1, v2};
  std::sort(_vertices.begin(), _vertices.end(),
            [](const Vertex *lhs, const Vertex *rhs) { return *lhs < *rhs; });
}

Triangle::Triangle(
    Vertex &vertexOne,
    Vertex &vertexTwo,
    Vertex &vertexThree)
    : _vertices({&vertexOne, &vertexTwo, &vertexThree})
{
  PRECICE_ASSERT(vertexOne.getDimensions() == vertexTwo.getDimensions(),
                 vertexOne.getDimensions(), vertexTwo.getDimensions());
  PRECICE_ASSERT(vertexTwo.getDimensions() == vertexThree.getDimensions(),
                 vertexTwo.getDimensions(), vertexThree.getDimensions());
  std::sort(_vertices.begin(), _vertices.end(),
            [](const Vertex *lhs, const Vertex *rhs) { return *lhs < *rhs; });
}

double Triangle::getArea() const
{
  return math::geometry::triangleArea(vertex(0).getCoords(), vertex(1).getCoords(), vertex(2).getCoords());
}

Eigen::VectorXd Triangle::computeNormal() const
{
  Eigen::Vector3d vectorA = (vertex(1).getCoords() - vertex(0).getCoords()) / 2.0;
  Eigen::Vector3d vectorB = (vertex(1).getCoords() - vertex(0).getCoords()) / 2.0;

  // Compute cross-product of vector A and vector B
  return vectorA.cross(vectorB).normalized();
}

int Triangle::getDimensions() const
{
  return _vertices[0]->getDimensions();
}

const Eigen::VectorXd Triangle::getCenter() const
{
  return (_vertices[0]->getCoords() + _vertices[1]->getCoords() + _vertices[2]->getCoords()) / 3.0;
}

double Triangle::getEnclosingRadius() const
{
  auto center = getCenter();
  return std::max({(center - _vertices[0]->getCoords()).norm(),
                   (center - _vertices[1]->getCoords()).norm(),
                   (center - _vertices[2]->getCoords()).norm()});
}

bool Triangle::operator==(const Triangle &other) const
{
  return std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin(),
                             [](const Vertex *e1, const Vertex *e2) { return *e1 == *e2; });
}

bool Triangle::operator!=(const Triangle &other) const
{
  return !(*this == other);
}

std::ostream &operator<<(std::ostream &os, const Triangle &t)
{
  using utils::eigenio::wkt;
  return os << "POLYGON (("
            << t.vertex(0).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(1).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(2).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(0).getCoords().transpose().format(wkt()) << "))";
}

} // namespace precice::mesh

#include "Triangle.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <boost/concept/assert.hpp>
#include <boost/range/concepts.hpp>
#include "math/differences.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Triangle::iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Triangle::const_iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<Triangle>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<const Triangle>) );

Triangle::Triangle(
    Edge &edgeOne,
    Edge &edgeTwo,
    Edge &edgeThree,
    int   id)
    : _edges({&edgeOne, &edgeTwo, &edgeThree}),
      _id(id),
      _normal(Eigen::VectorXd::Zero(edgeOne.getDimensions()))
{
  PRECICE_ASSERT(edgeOne.getDimensions() == edgeTwo.getDimensions(),
                 edgeOne.getDimensions(), edgeTwo.getDimensions());
  PRECICE_ASSERT(edgeTwo.getDimensions() == edgeThree.getDimensions(),
                 edgeTwo.getDimensions(), edgeThree.getDimensions());
  PRECICE_ASSERT(getDimensions() == 3, getDimensions());

  // Determine vertex map
  Vertex &v0 = edge(0).vertex(0);
  Vertex &v1 = edge(0).vertex(1);

  if (&edge(1).vertex(0) == &v0) {
    _vertexMap[0] = 1;
    _vertexMap[1] = 0;
  } else if (&edge(1).vertex(1) == &v0) {
    _vertexMap[0] = 1;
    _vertexMap[1] = 1;
  } else if (&edge(1).vertex(0) == &v1) {
    _vertexMap[0] = 0;
    _vertexMap[1] = 0;
  } else {
    PRECICE_ASSERT(&edge(1).vertex(1) == &v1);
    _vertexMap[0] = 0;
    _vertexMap[1] = 1;
  }

  if (_vertexMap[1] == 0) {
    if (&edge(2).vertex(0) == &edge(1).vertex(1)) {
      _vertexMap[2] = 0;
    } else {
      PRECICE_ASSERT(&edge(2).vertex(1) == &edge(1).vertex(1));
      _vertexMap[2] = 1;
    }
  } else if (_vertexMap[1] == 1) {
    if (&edge(2).vertex(0) == &edge(1).vertex(0)) {
      _vertexMap[2] = 0;
    } else {
      PRECICE_ASSERT(&edge(2).vertex(1) == &edge(1).vertex(0));
      _vertexMap[2] = 1;
    }
  }

  PRECICE_ASSERT(
      (&edge(0).vertex(_vertexMap[0]) != &edge(1).vertex(_vertexMap[1])) &&
          (&edge(0).vertex(_vertexMap[0]) != &edge(2).vertex(_vertexMap[2])) &&
          (&edge(1).vertex(_vertexMap[1]) != &edge(2).vertex(_vertexMap[2])),
      "Triangle vertices are not unique!");
}

double Triangle::getArea() const
{
  Eigen::Vector3d vectorA = edge(1).vertex(1).getCoords() - edge(1).vertex(0).getCoords();
  Eigen::Vector3d vectorB = edge(0).vertex(1).getCoords() - edge(0).vertex(0).getCoords();
  // Compute cross-product of vector A and vector B
  auto normal = vectorA.cross(vectorB);
  return (0.5 * normal.norm());
}

const Eigen::VectorXd Triangle::computeNormal(bool flip)
{
  Eigen::Vector3d vectorA = edge(1).getCenter() - edge(0).getCenter();
  Eigen::Vector3d vectorB = edge(2).getCenter() - edge(0).getCenter();
  // Compute cross-product of vector A and vector B
  auto normal = vectorA.cross(vectorB);
  if (flip) {
    normal *= -1.0; // Invert direction if counterclockwise
  }
  _normal = normal.normalized();
  return normal;
}

int Triangle::getDimensions() const
{
  return _edges[0]->getDimensions();
}

const Eigen::VectorXd &Triangle::getNormal() const
{
  return _normal;
}

const Eigen::VectorXd Triangle::getCenter() const
{
  return (_edges[0]->getCenter() + _edges[1]->getCenter() + _edges[2]->getCenter()) / 3.0;
}

double Triangle::getEnclosingRadius() const
{
  auto center = getCenter();
  return std::max({(center - vertex(0).getCoords()).norm(),
                   (center - vertex(1).getCoords()).norm(),
                   (center - vertex(2).getCoords()).norm()});
}

bool Triangle::operator==(const Triangle &other) const
{
  return math::equals(_normal, other._normal) &&
         std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin(),
                             [](const Edge *e1, const Edge *e2) { return *e1 == *e2; });
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

} // namespace mesh
} // namespace precice

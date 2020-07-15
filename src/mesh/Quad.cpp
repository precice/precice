#include "Quad.hpp"
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

BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Quad::iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Quad::const_iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<Quad>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<const Quad>) );

Quad::Quad(
    Edge &edgeOne,
    Edge &edgeTwo,
    Edge &edgeThree,
    Edge &edgeFour,
    int   id)
    : _edges({&edgeOne, &edgeTwo, &edgeThree, &edgeFour}),
      _id(id),
      _normal(Eigen::VectorXd::Zero(edgeOne.getDimensions()))
{
  PRECICE_ASSERT(edgeOne.getDimensions() == edgeTwo.getDimensions(),
                 edgeOne.getDimensions(), edgeTwo.getDimensions());
  PRECICE_ASSERT(edgeTwo.getDimensions() == edgeThree.getDimensions(),
                 edgeTwo.getDimensions(), edgeThree.getDimensions());
  PRECICE_ASSERT(edgeThree.getDimensions() == edgeFour.getDimensions(),
                 edgeThree.getDimensions(), edgeFour.getDimensions());
  PRECICE_ASSERT(getDimensions() == 3, getDimensions());

  // Determine vertex map
  Vertex &v0 = edge(0).vertex(0);
  Vertex &v1 = edge(0).vertex(1);

  // Check for edges 0 and 1 which vertex establishes connection
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

  // Check for edges 1 and 2 which vertex establishes connection
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

  // Check for edges 2 and 3 which vertex establishes connection
  if (_vertexMap[2] == 0) {
    if (&edge(3).vertex(0) == &edge(2).vertex(1)) {
      _vertexMap[3] = 0;
    } else {
      PRECICE_ASSERT(&edge(3).vertex(1) == &edge(2).vertex(1));
      _vertexMap[3] = 1;
    }
  } else if (_vertexMap[2] == 1) {
    if (&edge(3).vertex(0) == &edge(2).vertex(0)) {
      _vertexMap[3] = 0;
    } else {
      PRECICE_ASSERT(&edge(3).vertex(1) == &edge(2).vertex(0));
      _vertexMap[3] = 1;
    }
  }

  PRECICE_ASSERT(&vertex(0) != &vertex(1));
  PRECICE_ASSERT(&vertex(0) != &vertex(2));
  PRECICE_ASSERT(&vertex(0) != &vertex(3));
  PRECICE_ASSERT(&vertex(1) != &vertex(2));
  PRECICE_ASSERT(&vertex(1) != &vertex(3));
  PRECICE_ASSERT(&vertex(2) != &vertex(3));

  PRECICE_ASSERT((_vertexMap[0] == 0) || (_vertexMap[0] == 1), _vertexMap[0]);
  PRECICE_ASSERT((_vertexMap[1] == 0) || (_vertexMap[1] == 1), _vertexMap[1]);
  PRECICE_ASSERT((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[2]);
  PRECICE_ASSERT((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[3]);
}

const Eigen::VectorXd Quad::computeNormal(bool flip)
{
  // Two triangles are thought by splitting the quad from vertex 0 to 2.
  // The cross prodcut of the outer edges of the triangles is used to compute
  // the normal direction and area of the triangles. The direction must be
  // the same, while the areas differ in general. The normals are added up
  // and divided by 2 to get the area of the overall quad, since the length
  // does correspond to the parallelogram spanned by the vectors of the
  // cross product, which is twice the area of the corresponding triangles.
  Eigen::Vector3d vectorA = vertex(2).getCoords() - vertex(1).getCoords();
  Eigen::Vector3d vectorB = vertex(0).getCoords() - vertex(1).getCoords();
  // Compute cross-product of vector A and vector B
  auto normal = vectorA.cross(vectorB);

  vectorA               = vertex(0).getCoords() - vertex(3).getCoords();
  vectorB               = vertex(2).getCoords() - vertex(3).getCoords();
  auto normalSecondPart = vectorA.cross(vectorB);

  PRECICE_ASSERT(math::equals(normal.normalized(), normalSecondPart.normalized()),
                 normal, normalSecondPart);
  normal += normalSecondPart;
  normal *= 0.5;

  if (flip) {
    normal *= -1.0; // Invert direction if counterclockwise
  }
  _normal = normal.normalized();
  return normal;
}

int Quad::getDimensions() const
{
  return _edges[0]->getDimensions();
}

const Eigen::VectorXd &Quad::getNormal() const
{
  return _normal;
}

const Eigen::VectorXd Quad::getCenter() const
{
  return (_edges[0]->getCenter() + _edges[1]->getCenter() + _edges[2]->getCenter() + _edges[3]->getCenter()) / 4;
}

double Quad::getEnclosingRadius() const
{
  auto center = getCenter();
  return std::max({(center - vertex(0).getCoords()).norm(),
                   (center - vertex(1).getCoords()).norm(),
                   (center - vertex(2).getCoords()).norm(),
                   (center - vertex(3).getCoords()).norm()});
}

bool Quad::operator==(const Quad &other) const
{
  return math::equals(_normal, other._normal) &&
         std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin(),
                             [](const Edge *e1, const Edge *e2) { return *e1 == *e2; });
}

bool Quad::operator!=(const Quad &other) const
{
  return !(*this == other);
}

std::ostream &operator<<(std::ostream &os, const Quad &q)
{
  using utils::eigenio::wkt;
  return os << "POLYGON (("
            << q.vertex(0).getCoords().transpose().format(wkt()) << ", "
            << q.vertex(1).getCoords().transpose().format(wkt()) << ", "
            << q.vertex(2).getCoords().transpose().format(wkt()) << ", "
            << q.vertex(3).getCoords().transpose().format(wkt()) << ", "
            << q.vertex(0).getCoords().transpose().format(wkt()) << "))";
}

} // namespace mesh
} // namespace precice

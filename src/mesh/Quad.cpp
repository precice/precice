#include "Quad.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include <boost/range/concepts.hpp>

namespace precice
{
namespace mesh
{

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
  assertion(edgeOne.getDimensions() == edgeTwo.getDimensions(),
            edgeOne.getDimensions(), edgeTwo.getDimensions());
  assertion(edgeTwo.getDimensions() == edgeThree.getDimensions(),
            edgeTwo.getDimensions(), edgeThree.getDimensions());
  assertion(edgeThree.getDimensions() == edgeFour.getDimensions(),
            edgeThree.getDimensions(), edgeFour.getDimensions());
  assertion(getDimensions() == 3, getDimensions());

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
    assertion(&edge(1).vertex(1) == &v1);
    _vertexMap[0] = 0;
    _vertexMap[1] = 1;
  }

  // Check for edges 1 and 2 which vertex establishes connection
  if (_vertexMap[1] == 0) {
    if (&edge(2).vertex(0) == &edge(1).vertex(1)) {
      _vertexMap[2] = 0;
    } else {
      assertion(&edge(2).vertex(1) == &edge(1).vertex(1));
      _vertexMap[2] = 1;
    }
  } else if (_vertexMap[1] == 1) {
    if (&edge(2).vertex(0) == &edge(1).vertex(0)) {
      _vertexMap[2] = 0;
    } else {
      assertion(&edge(2).vertex(1) == &edge(1).vertex(0));
      _vertexMap[2] = 1;
    }
  }

  // Check for edges 2 and 3 which vertex establishes connection
  if (_vertexMap[2] == 0) {
    if (&edge(3).vertex(0) == &edge(2).vertex(1)) {
      _vertexMap[3] = 0;
    } else {
      assertion(&edge(3).vertex(1) == &edge(2).vertex(1));
      _vertexMap[3] = 1;
    }
  } else if (_vertexMap[2] == 1) {
    if (&edge(3).vertex(0) == &edge(2).vertex(0)) {
      _vertexMap[3] = 0;
    } else {
      assertion(&edge(3).vertex(1) == &edge(2).vertex(0));
      _vertexMap[3] = 1;
    }
  }

  assertion(&vertex(0) != &vertex(1));
  assertion(&vertex(0) != &vertex(2));
  assertion(&vertex(0) != &vertex(3));
  assertion(&vertex(1) != &vertex(2));
  assertion(&vertex(1) != &vertex(3));
  assertion(&vertex(2) != &vertex(3));

  assertion((_vertexMap[0] == 0) || (_vertexMap[0] == 1), _vertexMap[0]);
  assertion((_vertexMap[1] == 0) || (_vertexMap[1] == 1), _vertexMap[1]);
  assertion((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[2]);
  assertion((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[3]);
}

int Quad::getDimensions() const
{
  return _edges[0]->getDimensions();
}

void Quad::setEnclosingRadius(double radius)
{
  _enclosingRadius = radius;
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
  return _enclosingRadius;
}

bool Quad::operator==(const Quad& other) const
{
    return math::equals(_normal, other._normal) &&
        std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin(),
                [](const Edge* e1, const Edge* e2){return *e1 == *e2;});
}

bool Quad::operator!=(const Quad& other) const
{
  return !(*this == other);
}

std::ostream& operator<<(std::ostream& os, const Quad& q)
{
    os << "POLYGON ((";
    for (int i = 0; i < 4; i++){
        os << q.vertex(i).getCoords().transpose();
        if (i < 3)
            os << ", ";
    }
    return os <<", " << q.vertex(0).getCoords().transpose() << "))";
}

} // namespace mesh
} // namespace precice

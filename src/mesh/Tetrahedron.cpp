#include "Tetrahedron.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <boost/concept/assert.hpp>
#include <boost/range/concepts.hpp>
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Tetrahedron::iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessIteratorConcept<Tetrahedron::const_iterator>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<Tetrahedron>) );
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<const Tetrahedron>) );

Tetrahedron::Tetrahedron(
    Vertex &      vertexOne,
    Vertex &      vertexTwo,
    Vertex &      vertexThree,
    Vertex &      vertexFour,
    TetrahedronID id)
    : _vertices({&vertexOne, &vertexTwo, &vertexThree, &vertexFour}),
      _id(id)
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
  using utils::eigenio::wkt;
  return os << "POLYGON (("
            << t.vertex(0).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(1).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(2).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(3).getCoords().transpose().format(wkt()) << ", "
            << t.vertex(0).getCoords().transpose().format(wkt()) << "))";
}

} // namespace mesh
} // namespace precice

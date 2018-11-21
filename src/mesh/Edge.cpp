#include "Edge.hpp"

namespace precice {
namespace mesh {

Edge:: Edge
(
  Vertex& vertexOne,
  Vertex& vertexTwo,
  int     id )
:
  _vertices( {&vertexOne, &vertexTwo} ),
  _id ( id ),
  _normal ( Eigen::VectorXd::Constant(vertexOne.getDimensions(), 0.0) )
{
  assertion ( vertexOne.getDimensions() == vertexTwo.getDimensions(),
              vertexOne.getDimensions(), vertexTwo.getDimensions() );
}

int Edge:: getID () const
{
  return _id;
}

const Eigen::VectorXd Edge::getCenter () const
{
  return 0.5 * (_vertices[0]->getCoords() + _vertices[1]->getCoords());
}

double Edge:: getEnclosingRadius () const
{
  return (_vertices[0]->getCoords() - getCenter()).norm();
}

bool Edge::operator==(const Edge& other) const
{
    return math::equals(_normal, other._normal) &&
        std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin(),
                [](const Vertex* a, const Vertex* b){return *a == *b;});
}
bool Edge::operator!=(const Edge& other) const
{
  return !(*this == other);
}

std::ostream& operator<<(std::ostream& stream, const Edge& edge){
    stream << "LINESTRING (";
    for (int i = 0; i < 2; i++){
        stream << edge.vertex(i).getCoords().transpose();
        if (i < 1)
            stream << ", ";
    }
    return stream << ")";
}

}} // namespace precice, mesh

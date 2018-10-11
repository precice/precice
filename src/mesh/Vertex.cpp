#include "Vertex.hpp"

namespace precice {
namespace mesh {

int Vertex:: getDimensions() const
{
  return _coords.size();
}

const Eigen::VectorXd& Vertex::getNormal () const
{
  return _normal;
}

int Vertex:: getGlobalIndex() const {
  return _globalIndex;
}

void Vertex:: setGlobalIndex(int globalIndex){
  _globalIndex = globalIndex;
}

bool Vertex:: isOwner() const {
  return _owner;
}

void Vertex:: setOwner(bool owner){
  _owner = owner;
}

bool Vertex:: isTagged() const {
  return _tagged;
}

void Vertex:: tag() {
  _tagged = true;
}


std::ostream & operator<<(std::ostream &os, Vertex const & v)
{
  // transpose, so output is on one line
  return os << "Vertex " << v.getID() << ": " << v.getCoords().transpose() << "\n"; 
}

}} // namespace precice, mesh

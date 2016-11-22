#include "Vertex.hpp"
#include "utils/ManageUniqueIDs.hpp"

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

const Mesh* Vertex:: mesh () const
{
  return _mesh;
}

Mesh* Vertex:: mesh ()
{
  return _mesh;
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

}} // namespace precice, mesh

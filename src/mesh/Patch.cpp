#include "Patch.hpp"
#include "Vertex.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

//void Patch::setPatchID(int patchID)
//{
//  _patchid = patchID;
//}

Patch::Patch(
    const std::string &name,
    int                id)
    : _name(name),
      _id(id)
{

}

Patch::patchVertexContainer &Patch::patchVertices()
{
  return _vertices;
}

const Patch::patchVertexContainer &Patch::patchVertices() const
{
  return _vertices;
}

const std::string &Patch::getName() const
{
  return _name;
}

int Patch::getID() const
{
  return _id;
}



} // namespace mesh
} // namespace precice

#include "Vertex.hpp"
#include "utils/EigenIO.hpp"

namespace precice {
namespace mesh {

int Vertex::getDimensions() const
{
  return _coords.size();
}

const Eigen::VectorXd &Vertex::getNormal() const
{
  return _normal;
}

int Vertex::getGlobalIndex() const
{
  return _globalIndex;
}

void Vertex::setGlobalIndex(int globalIndex)
{
  _globalIndex = globalIndex;
}

bool Vertex::isOwner() const
{
  return _owner;
}

void Vertex::setOwner(bool owner)
{
  _owner = owner;
}

bool Vertex::isTagged() const
{
  return _tagged;
}

void Vertex::tag()
{
  _tagged = true;
}

void Vertex::setPatchID(int patchid)
{
  _patchid = patchid;
}

void Vertex::setPatchName(std::string patchname)
{
  _patchname = patchname;
}

std::ostream &operator<<(std::ostream &os, Vertex const &v)
{
  return os << "POINT (" << v.getCoords().transpose().format(utils::eigenio::wkt()) << ')';
}

} // namespace mesh
} // namespace precice

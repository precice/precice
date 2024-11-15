#include "Vertex.hpp"
#include <Eigen/Core>
#include "utils/EigenIO.hpp"

namespace precice::mesh {

int Vertex::getDimensions() const
{
  return _dim;
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

std::ostream &operator<<(std::ostream &os, Vertex const &v)
{
  return os << "POINT (" << v.getCoords().transpose().format(utils::eigenio::wkt()) << ')';
}

} // namespace precice::mesh

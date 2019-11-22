#include "Merge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace mesh {

Group& Merge:: content()
{
  return _merged;
}

const Group& Merge:: content() const
{
  return _merged;
}

}} // namespace precice, mesh


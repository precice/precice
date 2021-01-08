#pragma once

#include "mesh/BoundingBox.hpp"
#include "query/RTree.hpp"

namespace precice {
namespace mesh {

inline query::RTreeBox toRTreeBox(BoundingBox const &bb)
{
  return query::RTreeBox{bb.minCorner(), bb.maxCorner()};
}

} // namespace mesh
} // namespace precice

#pragma once

#include "mesh/BoundingBox.hpp"
#include "mesh/RTree.hpp"

namespace precice {
namespace mesh {

inline RTreeBox toRTreeBox(BoundingBox const &bb)
{
  return RTreeBox{bb.minCorner(), bb.maxCorner()};
}

} // namespace mesh
} // namespace precice

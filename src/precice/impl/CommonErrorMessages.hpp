#pragma once

#include <string>

namespace precice {

namespace impl {

inline std::string errorInvalidVertexID(int vid)
{
  return "The given VertexID \"" + std::to_string(vid) + "\" is invalid. Check that it originated from a call to setMeshVertex() or setMeshVertices().";
}

inline std::string errorInvalidEdgeID(int eid)
{
  return "The given EdgeID \"" + std::to_string(eid) + "\" is invalid. Check that it originated from a call to setMeshEdge().";
}

static constexpr auto errorInvalidVertexIDRange = "The given range of VertexIDs contains invalid IDs at offsets [{},{}]. Check that they originated from calls to setMeshVertex() or setMeshVertices().";

} // namespace impl

} // namespace precice

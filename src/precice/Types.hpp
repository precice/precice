#pragma once

namespace precice {

/**
 * Type used identify a vertex on a mesh.
 *
 * These IDs are always bound to a mesh. Using a VertexID on another mesh is not allowed and not guaranteed to result in an error.
 * The only valid way of receiving these IDs are by using \ref setMeshVertex(), \ref setMeshVertices(), \ref getMeshVertexIDsAndCoordinates().
 *
 * There is no guarantee on the values of vertex IDs. **Do not** rely on them to start at zero or be continuous.
 */
using VertexID = int;

} // namespace precice

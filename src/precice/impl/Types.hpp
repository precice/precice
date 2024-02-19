#pragma once

#include "precice/Types.hpp"

namespace precice {

/**
 * Type used for the IDs of edges
 */
using EdgeID = int;

/**
 * Type used for the IDs of triangles
 */
using TriangleID = int;

/**
 * Type used for the IDs of tetrahedra
 */
using TetrahedronID = int;

/**
 * Type used for the IDs of data (incl. gradient data)
 */
using DataID = int;

/**
 * Type used for the IDs of meshes
 */
using MeshID = int;

/**
 * Type used to represent local ranks of a parallel participant.
 *
 * The primary participant has Rank 0.
 */
using Rank = int;

} // namespace precice

#pragma once

#include <Eigen/Core>
#include <iostream>

#include "math/differences.hpp"

namespace precice {
namespace mesh {

/// Patch of a vertex for a mesh.
class Patch {
public:
  using patchVertexContainer   = std::deque<Vertex>;

  /**
   * @brief Constructor.
   *
   * @param[in] name Unique name of the mesh.
   * @param[in] dimensions Dimensionalty of the mesh.
   * @param[in] flipNormals Inverts the standard direction of normals.
   * @param[in] id The id of this mesh
   */
  Patch(
      const std::string &name,
      int                id);
      
  /// Returns modifieable container holding all vertices.
  patchVertexContainer &patchVertices();

  /// Returns const container holding all vertices.
  const patchVertexContainer &patchVertices() const;

  /// Returns the name of the mesh, as set in the config file.
  const std::string &getName() const;

  /// Returns the base ID of the mesh.
  int getID() const;


private:

  /// Patch id that a vertex belongs to
  //int _patchid;

  /// Coordinates of the vertex.
  //Eigen::VectorXd _coords;

  /// Normal of the vertex.
  //Eigen::VectorXd _normal;

  /// global (unique) index for parallel simulations
  //int _globalIndex = -1;

  /// true if this processors is the owner of the vertex (for parallel simulations)
  //bool _owner = true;

  /// true if this vertex is tagged for partition
  //bool _tagged = false;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION
/*template <typename VECTOR_T>
Patch::Patch(
    int             id,
    int             patchid)
    : _id(id),
      _patchid(patchid)
{
}


inline int Patch::getPatchID() const
{
  return _patchid;

*/

} // namespace mesh
} // namespace precice

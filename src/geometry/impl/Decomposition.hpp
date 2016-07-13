#pragma once

#include "tarch/logging/Log.h"
#include "mapping/SharedPointer.hpp"
#include "utils/MasterSlave.hpp"
#include <map>
#include <vector>
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {
namespace impl {

/**
 * @brief Abstract base class for all decompositions.
 *
 * A decompositions decomposes a geometry among distributed memory.
 * This can be a received CommnuicatedGeometry, but also a built-in or an imported
 * Geometry. Before decomposition the complete underlying mesh is held by the master process.
 * After decomposition, every processor, slaves and master, hold a local version of the mesh.
 */
class Decomposition
{
public:

  /**
   * @brief Constructor.
   *
   * @param dimensions: spatial dimension of the geometry
   *
   */
  Decomposition ( int dimensions, double safetyFactor );

  /**
   * @brief Destructor.
   */
  virtual ~Decomposition() {}

  /**
   * @brief Decomposes the geometry.
   */
  virtual void decompose ( mesh::Mesh& seed ) =0;

  /**
   * @brief Sets an offset from zero for the geometry.
   */

  void setBoundingFromMapping(mapping::PtrMapping mapping);

  void setBoundingToMapping(mapping::PtrMapping mapping);

protected:

  void computeBoundingMappings();

  void clearBoundingMappings();

  std::vector<int> filterMesh(mesh::Mesh& seed, mesh::Mesh& filteredMesh);

  /// Returns true if a vertex contributes. If false, the vertex can be erased.
  bool doesVertexContribute(const mesh::Vertex& vertex);

  void mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb);

  /**
   *  @brief spatial dimension of the geometry
   */
  int _dimensions;

  mapping::PtrMapping _boundingFromMapping;

  mapping::PtrMapping _boundingToMapping;

  mesh::Mesh::BoundingBox _bb;

  double _safetyGap;

  double _safetyFactor;

  bool _filterByMapping;


private:

  static tarch::logging::Log _log;

};

}}} // namespace precice, geometry, impl


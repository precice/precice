// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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
  Decomposition ( int dimensions);

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

  /**
   * @brief send feedback about decomposition back to master
   */
  void feedback(
      mesh::Mesh& seed,
      std::map<int,std::vector<int> >& boundingVertexDistribution,
      std::vector<int>& filteredVertexPositions);

  /// Returns true if a vertex contributes. If false, the vertex can be erased.
  virtual bool doesVertexContribute(const mesh::Vertex& vertex)=0;

  /**
   *  @brief spatial dimension of the geometry
   */
  int _dimensions;

  mapping::PtrMapping _boundingFromMapping;

  mapping::PtrMapping _boundingToMapping;


private:

  static tarch::logging::Log _log;

};

}}} // namespace precice, geometry, impl


#pragma once

#include "Decomposition.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/MasterSlave.hpp"
#include "tarch/logging/Log.h"
#include <map>
#include <vector>

namespace precice {
namespace geometry {
namespace impl {

/**
 * @brief Decomposes a geometry resp. mesh by broadcasting the complete mesh
 * and a post-filter, afterwards, on each slave.
 */
class BroadcastFilterDecomposition : public Decomposition
{
public:

  BroadcastFilterDecomposition (
    int    dimensions, double safetyFactor );

  virtual ~BroadcastFilterDecomposition() {}

private:

  /**
   * @brief Decomposes the geometry.
   */
  void decompose(
    mesh::Mesh& seed);

  void broadcast(
    mesh::Mesh& seed);

  void filter(
    mesh::Mesh& seed,
    std::vector<int>& filteredVertexPositions);

  /**
   * @brief send feedback about decomposition back to master
   */
  void feedback(
    mesh::Mesh& seed,
    std::vector<int>& filteredVertexPositions);

  /// Logging device.
  static tarch::logging::Log _log;

};

}}} // namespace precice, geometry, filter

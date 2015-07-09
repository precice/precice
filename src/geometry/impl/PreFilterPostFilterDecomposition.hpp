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
 * @brief Decomposes a geometry resp. mesh by applying a pre-filter first (bounding box) on the master
 * and a post-filter, afterwards, on each slave.
 */
class PreFilterPostFilterDecomposition : public Decomposition
{
public:

  PreFilterPostFilterDecomposition (
    int    dimensions,
    double safetyFactor);

  virtual ~PreFilterPostFilterDecomposition() {}

private:

  /**
   * @brief Decomposes the geometry.
   */
  void decompose(
    mesh::Mesh& seed);

  void preFilter(
    mesh::Mesh& seed,
    std::map<int,std::vector<int> >& boundingVertexDistribution);

  void postFilter(
    mesh::Mesh& seed,
    std::map<int,std::vector<int> >& boundingVertexDistribution,
    std::vector<int>& filteredVertexPositions);


  void mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb);

  /// Returns true if a vertex contributes. If false, the vertex can be erased.
  bool doesVertexContribute(const mesh::Vertex& vertex);

  /// Logging device.
  static tarch::logging::Log _log;

  mesh::Mesh::BoundingBox _bb;

  double _safetyGap;

  double _safetyFactor;

  bool _filterByMapping;
};

}}} // namespace precice, geometry, filter

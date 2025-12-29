#pragma once

#include <vector>
#include "MappingContext.hpp"
#include "SharedPointer.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include "partition/SharedPointer.hpp"

namespace precice::impl {

/// Stores a mesh and related objects and data.
struct MeshContext {
  /** Upgrades the mesh requirement to a more specific level.
   * @param[in] requirement The requirement to upgrade to.
   */
  void require(mapping::Mapping::MeshRequirement requirement);

  /// Mesh holding the geometry data structure.
  mesh::PtrMesh mesh;

  /// Determines which mesh type has to be provided by the accessor.
  mapping::Mapping::MeshRequirement meshRequirement = mapping::Mapping::MeshRequirement::UNDEFINED;

  /// Name of participant that creates the mesh.
  std::string receiveMeshFrom;

  /// bounding box to speed up decomposition of received mesh is increased by this safety factor
  double safetyFactor = -1;

  /// In case a mapping done by the solver is favored over a preCICE mapping, apply user-defined
  /// bounding-boxes.
  bool allowDirectAccess = false;

  /// setMeshAccessRegion may be called multiple times per mesh(context).
  /// putting this into the mesh context means that we can only call
  /// this once, regardless of combinations with just-in-time mappings
  /// or multiple such mappings. If multiples are desired, we would
  /// need to shift this into the MappingContext
  std::vector<mesh::BoundingBox> userDefinedAccessRegions;

  /// True, if accessor does create the mesh.
  bool provideMesh = false;

  /// type of geometric filter
  partition::ReceivedPartition::GeometricFilter geoFilter = partition::ReceivedPartition::GeometricFilter::UNDEFINED;

  /// Partition creating the parallel decomposition of the mesh
  partition::PtrPartition partition;

  /// Mappings used when mapping data from the mesh. Can be empty.
  std::vector<MappingContext> fromMappingContexts;

  /// Mappings used when mapping data to the mesh. Can be empty.
  std::vector<MappingContext> toMappingContexts;

  /// Checks, that all vertices are within the user-defined access region and throws
  /// an error if vertices are not. The function does not return the result to be able
  /// to log the actual outliers
  void checkVerticesInsideAccessRegion(precice::span<const double> coordinates, const int meshDim, std::string_view functionName) const;

  std::vector<std::reference_wrapper<const mesh::Vertex>> filterVerticesToLocalAccessRegion(bool requiresBB) const;

  void clearMappings()
  {
    for (auto &mc : fromMappingContexts) {
      mc.mapping->clear();
    }
    for (auto &mc : toMappingContexts) {
      mc.mapping->clear();
    }
  }

private:
  mutable logging::Logger _log{"impl::MeshContext"};
};

inline void MeshContext::require(mapping::Mapping::MeshRequirement requirement)
{
  meshRequirement = std::max(meshRequirement, requirement);
}

inline void MeshContext::checkVerticesInsideAccessRegion(precice::span<const double> coordinates, const int meshDim, std::string_view functionName) const
{

  if (!userDefinedAccessRegions.empty()) {
    const auto                        nVertices = (coordinates.size() / meshDim);
    Eigen::Map<const Eigen::MatrixXd> C(coordinates.data(), meshDim, nVertices);
    Eigen::VectorXd                   minCoeffs = C.rowwise().minCoeff();
    Eigen::VectorXd                   maxCoeffs = C.rowwise().maxCoeff();

    const auto vertexInside = [&](const mesh::BoundingBox &bb) {
      return (minCoeffs.array() >= bb.minCorner().array()).all() && (maxCoeffs.array() <= bb.maxCorner().array()).all();
    };
    const bool insideAnyRegion = std::any_of(userDefinedAccessRegions.begin(),
                                             userDefinedAccessRegions.end(),
                                             vertexInside);

    PRECICE_CHECK(insideAnyRegion, "The provided coordinates in \"{}\" are not within the access region defined with \"setMeshAccessRegion()\". "
                                   "Minimum corner of the provided values is (x,y,z) = ({}). "
                                   "Maximum corner of the provided values is (x,y,z) = ({}). ",
                  functionName, minCoeffs, maxCoeffs);
    C.colwise().maxCoeff();
  }
}

inline std::vector<std::reference_wrapper<const mesh::Vertex>> MeshContext::filterVerticesToLocalAccessRegion(bool requiresBB) const
{
  std::vector<std::reference_wrapper<const mesh::Vertex>> filteredVertices;
  const bool                                              regionDefined = !userDefinedAccessRegions.empty();

  const auto insideOneRegion = [&](const mesh::Vertex &v) {
    return std::any_of(userDefinedAccessRegions.begin(),
                       userDefinedAccessRegions.end(),
                       [&](const auto &bb) { return bb.contains(v); });
  };

  for (const auto &v : mesh->vertices()) {
    if (regionDefined ? insideOneRegion(v) // at least one box contains the vertex
                      : !requiresBB)       // no region list; keep everything
    {
      filteredVertices.emplace_back(v);
    }
  }
  return filteredVertices;
}

} // namespace precice::impl

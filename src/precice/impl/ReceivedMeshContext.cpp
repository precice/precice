#include "ReceivedMeshContext.hpp"

#include <Eigen/Core>

#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "utils/assertion.hpp"

namespace precice::impl {

void ReceivedMeshContext::checkVerticesInsideAccessRegion(precice::span<const double> coordinates, int meshDim, std::string_view functionName) const
{
  if (!userDefinedAccessRegion) {
    return;
  }
  const auto                        nVertices = (coordinates.size() / meshDim);
  Eigen::Map<const Eigen::MatrixXd> C(coordinates.data(), meshDim, nVertices);
  Eigen::VectorXd                   minCoeffs = C.rowwise().minCoeff();
  Eigen::VectorXd                   maxCoeffs = C.rowwise().maxCoeff();
  bool                              minCheck  = (minCoeffs.array() >= userDefinedAccessRegion->minCorner().array()).all();
  bool                              maxCheck  = (maxCoeffs.array() <= userDefinedAccessRegion->maxCorner().array()).all();
  PRECICE_CHECK(minCheck && maxCheck, "The provided coordinates in \"{}\" are not within the access region defined with \"setMeshAccessRegion()\". "
                                      "Minimum corner of the provided values is (x,y,z) = ({}), the minimum corner of the access region box is (x,y,z) = ({}). "
                                      "Maximum corner of the provided values is (x,y,z) = ({}), the maximum corner of the access region box is (x,y,z) = ({}). ",
                functionName, minCoeffs, userDefinedAccessRegion->minCorner(), maxCoeffs, userDefinedAccessRegion->maxCorner());
  C.colwise().maxCoeff();
}

std::vector<std::reference_wrapper<const mesh::Vertex>> ReceivedMeshContext::filterVerticesToLocalAccessRegion(bool requiresBB) const
{
  std::vector<std::reference_wrapper<const mesh::Vertex>> filteredVertices;
  for (const auto &v : mesh->vertices()) {
    // either the vertex lies within the region OR the user-defined region is not strictly necessary
    if (userDefinedAccessRegion) {
      // region is defined: only add if the vertex is inside the region
      if (userDefinedAccessRegion->contains(v)) {
        filteredVertices.push_back(std::cref(v));
      }
    } else if (!requiresBB) {
      // region is not defined, so if filtering isn't required, add all vertices
      filteredVertices.push_back(std::cref(v));
    }
  }
  return filteredVertices;
}

} // namespace precice::impl

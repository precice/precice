#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <utility>

#include "logging/LogMacros.hpp"
#include "mapping/BarycentricBaseMapping.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/Polation.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "profiling/Event.hpp"
#include "query/Index.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

BarycentricBaseMapping::BarycentricBaseMapping(Constraint constraint, int dimensions)
    : Mapping(constraint, dimensions, false, Mapping::InitialGuessRequirement::None)
{
}

void BarycentricBaseMapping::clear()
{
  PRECICE_TRACE();
  _operations.clear();
  _hasComputedMapping = false;
}

namespace {

template <int n>
void mapTemplatedConsistent(const std::vector<Operation> &ops, const Eigen::VectorXd &in, Eigen::VectorXd &out)
{
  static_assert(n > 0 && n <= 3);
  // For each output vertex, compute the linear combination of input vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  for (const auto &op : ops) {
    if constexpr (n == 1) {
      // We use a direct index-loop when no dimension offset is required
      out[op.out] += op.weight * in[op.in];
    } else {
      // Segments of templated size are the fastest option when the data has multiple components
      out.segment<n>(op.out * n).noalias() += op.weight * in.segment<n>(op.in * n);
    }
  }
}

template <int n>
void mapTemplatedConservative(const std::vector<Operation> &ops, const Eigen::VectorXd &in, Eigen::VectorXd &out)
{
  static_assert(n > 0 && n <= 3);
  // For each input vertex, distribute the conserved data among the relevant output vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  for (const auto &op : ops) {
    if constexpr (n == 1) {
      // We use a direct index-loop when no dimension offset is required
      out[op.in] += op.weight * in[op.out];
    } else {
      // Segments of templated size are the fastest option when the data has multiple components
      out.segment<n>(op.in * n).noalias() += op.weight * in.segment<n>(op.out * n);
    }
  }
}

} // namespace

void BarycentricBaseMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.bbm.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_DEBUG("Map conservative using {}", getName());
  const int              dimensions = inData.dataDims;
  const Eigen::VectorXd &inValues   = inData.values;
  Eigen::VectorXd       &outValues  = outData;

  switch (dimensions) {
  case 1:
    mapTemplatedConservative<1>(_operations, inValues, outValues);
    return;
  case 2:
    mapTemplatedConservative<2>(_operations, inValues, outValues);
    return;
  case 3:
    mapTemplatedConservative<3>(_operations, inValues, outValues);
    return;
  default:
    PRECICE_UNREACHABLE("Implement for unknown dimension");
  }
}

void BarycentricBaseMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.bbm.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map {} using {}", (hasConstraint(CONSISTENT) ? "consistent" : "scaled-consistent"), getName());

  const int              dimensions = inData.dataDims;
  const Eigen::VectorXd &inValues   = inData.values;
  Eigen::VectorXd       &outValues  = outData;

  switch (dimensions) {
  case 1:
    mapTemplatedConsistent<1>(_operations, inValues, outValues);
    return;
  case 2:
    mapTemplatedConsistent<2>(_operations, inValues, outValues);
    return;
  case 3:
    mapTemplatedConsistent<3>(_operations, inValues, outValues);
    return;
  default:
    PRECICE_UNREACHABLE("Implement for unknown dimension");
  }
}

void BarycentricBaseMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.bbm.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Compute Mapping for Tagging");

  computeMapping();
  PRECICE_DEBUG("Tagging First Round");

  // Determine the Mesh to Tag
  mesh::PtrMesh origins;
  if (hasConstraint(CONSERVATIVE)) {
    origins = output();
  } else {
    origins = input();
  }

  // Gather all vertices to be tagged in a first phase.
  std::vector<bool> tagged(origins->nVertices(), false);
  for (const auto &op : _operations) {
    PRECICE_ASSERT(!math::equals(op.weight, 0.0));
    tagged[op.in] = true;
  }

  // Now tag all vertices to be tagged in the second phase.
  size_t ntagged = 0;
  for (size_t i = 0; i < origins->nVertices(); ++i) {
    if (tagged[i]) {
      origins->vertex(i).tag();
      ++ntagged;
    }
  }

  PRECICE_DEBUG("First Round Tagged {}/{} Vertices", ntagged, origins->nVertices());

  clear();
}

void BarycentricBaseMapping::addPolation(VertexID out, const Polation &p)
{
  for (const auto &we : p.getWeightedElements()) {
    if (!precice::math::equals(we.weight, 0.0)) {
      _operations.push_back({out, we.vertexID, we.weight});
    }
  }
}

void BarycentricBaseMapping::postProcessOperations()
{
  if (hasConstraint(CONSERVATIVE)) {
    // We order the operations to make them write sequentially
    std::sort(_operations.begin(), _operations.end(), [](const auto &lhs, const auto &rhs) {
      // weak order for a pair
      if (lhs.in < rhs.in) {
        return true;
      }
      if (rhs.in < lhs.in) {
        return false;
      }
      // lhs.in == rhs.in
      return lhs.out < rhs.out;
    });
  }
}

void BarycentricBaseMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NP mapping no operation needed here
}

} // namespace precice::mapping

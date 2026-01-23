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

void BarycentricBaseMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.bbm.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_DEBUG("Map conservative using {}", getName());
  const int              dimensions = inData.dataDims;
  const Eigen::VectorXd &inValues   = inData.values;
  Eigen::VectorXd       &outValues  = outData;

  // For each input vertex, distribute the conserved data among the relevant output vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  if (dimensions == 1) {
    // Use non-strided access for 1D data
    for (const auto &op : _operations) {
      outValues[op.in] += op.weight * inValues[op.out];
    }
  } else {
    for (const auto &op : _operations) {
      outValues.segment(op.in * dimensions, dimensions) += op.weight * inValues.segment(op.out * dimensions, dimensions);
    }
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

  // For each output vertex, compute the linear combination of input vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  if (dimensions == 1) {
    // Use non-strided access for 1D data
    for (const auto &op : _operations) {
      outValues[op.out] += op.weight * inValues[op.in];
    }
  } else {
    for (const auto &op : _operations) {
      outValues.segment(op.out * dimensions, dimensions) += op.weight * inValues.segment(op.in * dimensions, dimensions);
    }
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
  // max_count is used to shortcut if all vertices have been tagged.
  std::unordered_set<int> tagged;
  const std::size_t       max_count = origins->nVertices();

  for (const auto &op : _operations) {
    PRECICE_ASSERT(!math::equals(op.weight, 0.0));
    tagged.insert(op.in);
    // Shortcut if all vertices are tagged
    if (tagged.size() == max_count) {
      break;
    }
  }

  // Now tag all vertices to be tagged in the second phase.
  for (auto &v : origins->vertices()) {
    if (tagged.count(v.getID()) == 1) {
      v.tag();
    }
  }
  PRECICE_DEBUG("First Round Tagged {}/{} Vertices", tagged.size(), max_count);

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

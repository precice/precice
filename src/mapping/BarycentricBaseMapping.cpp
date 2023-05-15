#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <utility>

#include "logging/LogMacros.hpp"
#include "mapping/BarycentricBaseMapping.hpp"
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
    : Mapping(constraint, dimensions)
{
}

void BarycentricBaseMapping::clear()
{
  PRECICE_TRACE();
  _interpolations.clear();
  _hasComputedMapping = false;
}

void BarycentricBaseMapping::mapConservative(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);
  precice::profiling::Event e("map.bbm.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_DEBUG("Map conservative");
  PRECICE_ASSERT(_interpolations.size() == input()->vertices().size(),
                 _interpolations.size(), input()->vertices().size());
  const int              dimensions = input()->data(inputDataID)->getDimensions();
  const Eigen::VectorXd &inValues   = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outValues  = output()->data(outputDataID)->values();

  // For each input vertex, distribute the conserved data among the relevant output vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  for (size_t i = 0; i < input()->vertices().size(); i++) {
    const size_t inOffset = i * dimensions;
    const auto & elems    = _interpolations[i].getWeightedElements();
    for (const auto &elem : elems) {
      size_t outOffset = static_cast<size_t>(elem.vertexID) * dimensions;
      for (int dim = 0; dim < dimensions; dim++) {
        PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
        PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
        outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
      }
    }
  }
}

void BarycentricBaseMapping::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);
  precice::profiling::Event e("map.bbm.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map consistent");
  PRECICE_ASSERT(_interpolations.size() == output()->vertices().size(),
                 _interpolations.size(), output()->vertices().size());

  const int              dimensions = input()->data(inputDataID)->getDimensions();
  const Eigen::VectorXd &inValues   = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outValues  = output()->data(outputDataID)->values();

  // For each output vertex, compute the linear combination of input vertices
  // Do it for all dimensions (i.e. components if data is a vector)
  for (size_t i = 0; i < output()->vertices().size(); i++) {
    const auto &elems     = _interpolations[i].getWeightedElements();
    size_t      outOffset = i * dimensions;
    for (const auto &elem : elems) {
      const size_t inOffset = static_cast<size_t>(elem.vertexID) * dimensions;
      for (int dim = 0; dim < dimensions; dim++) {
        PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
        PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
        outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
      }
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
  const std::size_t       max_count = origins->vertices().size();

  for (const Polation &interpolation : _interpolations) {
    for (const auto &elem : interpolation.getWeightedElements()) {
      if (!math::equals(elem.weight, 0.0)) {
        tagged.insert(elem.vertexID);
      }
    }
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

void BarycentricBaseMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NP mapping no operation needed here
}

} // namespace precice::mapping

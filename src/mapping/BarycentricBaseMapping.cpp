
#include <Eigen/Core>
#include <algorithm>
#include <deque>
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
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

BarycentricBaseMapping::BarycentricBaseMapping(Constraint constraint, int dimensions)
    : Mapping(constraint, dimensions)
{
}

bool BarycentricBaseMapping::hasComputedMapping() const
{
  return _hasComputedMapping;
}

void BarycentricBaseMapping::clear()
{
  PRECICE_TRACE();
  _interpolations.clear();
  _hasComputedMapping = false;
}

void BarycentricBaseMapping::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map.bm.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  mesh::PtrData          inData    = input()->data(inputDataID);
  mesh::PtrData          outData   = output()->data(outputDataID);
  const Eigen::VectorXd &inValues  = inData->values();
  Eigen::VectorXd &      outValues = outData->values();
const  int dimensions = inData->getDimensions();
  PRECICE_ASSERT(dimensions == outData->getDimensions());

  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    PRECICE_DEBUG("Map conservative");
    PRECICE_ASSERT(_interpolations.size() == input()->vertices().size(),
                   _interpolations.size(), input()->vertices().size());
    for (size_t i = 0; i < input()->vertices().size(); i++) {
    const  size_t      inOffset = i * dimensions;
      const auto &elems    = _interpolations[i].getWeightedElements();
      for (const auto &elem : elems) {
        size_t outOffset = static_cast<size_t>(elem.vertexID) * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
          PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  } else {
    PRECICE_DEBUG("Map consistent");
    PRECICE_ASSERT(_interpolations.size() == output()->vertices().size(),
                   _interpolations.size(), output()->vertices().size());
    for (size_t i = 0; i < output()->vertices().size(); i++) {
      const auto &elems     = _interpolations[i].getWeightedElements();
      size_t      outOffset = i * dimensions;
      for (const auto &elem : elems) {
        size_t inOffset = static_cast<size_t>(elem.vertexID) * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
          PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
    if (hasConstraint(SCALEDCONSISTENT)) {
      scaleConsistentMapping(inputDataID, outputDataID);
    }
  }
}

void BarycentricBaseMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map.bbm.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
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

} // namespace mapping
} // namespace precice

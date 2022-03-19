#include "NearestProjectionMapping.hpp"

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

BarycentricBaseMapping::BarycentricBaseMapping(Constraint constraint, int dimensions) : Mapping(constraint, dimensions) {

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

  precice::utils::Event e("map.np.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  mesh::PtrData          inData    = input()->data(inputDataID);
  mesh::PtrData          outData   = output()->data(outputDataID);
  const Eigen::VectorXd &inValues  = inData->values();
  Eigen::VectorXd &      outValues = outData->values();
  //assign(outValues) = 0.0;
  int dimensions = inData->getDimensions();
  PRECICE_ASSERT(dimensions == outData->getDimensions());

  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    PRECICE_DEBUG("Map conservative");
    PRECICE_ASSERT(_interpolations.size() == input()->vertices().size(),
                   _interpolations.size(), input()->vertices().size());
    for (size_t i = 0; i < input()->vertices().size(); i++) {
      size_t      inOffset = i * dimensions;
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


} // namespace mapping
} // namespace precice

#include "NearestProjectionMapping.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <utility>

#include <boost/version.hpp>
#if BOOST_VERSION < 106600
#include <boost/function_output_iterator.hpp>
#else
#include <boost/iterator/function_output_iterator.hpp>
#endif

#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "query/Interpolation.hpp"
#include "query/RTree.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace precice {
extern bool syncMode;

namespace mapping {

NearestProjectionMapping::NearestProjectionMapping(
    Constraint constraint,
    int        dimensions)
    : Mapping(constraint, dimensions)
{
  if (constraint == CONSISTENT) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  } else {
    PRECICE_ASSERT(constraint == CONSERVATIVE, constraint);
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }
}

void NearestProjectionMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size(), output()->vertices().size());
  const std::string     baseEvent = "map.np.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  // Setup Direction of Mapping
  mesh::PtrMesh origins, searchSpace;
  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    origins     = output();
    searchSpace = input();
  } else {
    PRECICE_DEBUG("Compute conservative mapping");
    origins     = input();
    searchSpace = output();
  }

  const auto &fVertices = origins->vertices();
  const auto &tVertices = searchSpace->vertices();
  const auto &tEdges    = searchSpace->edges();

  _weights.resize(fVertices.size());

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  if (getDimensions() == 2) {
    if (!fVertices.empty() && tEdges.empty()) {
      PRECICE_WARN("2D Mesh \"" << searchSpace->getName() << "\" does not contain edges. Nearest projection mapping falls back to nearest neighbor mapping.");
    }

    // Lazy evaluation of the vertex index.
    // This is not necessary in the case of matching meshes.
    utils::statistics::DistanceAccumulator distanceStatistics;

    for (size_t i = 0; i < fVertices.size(); i++) {
      // Search for the origin inside the destination meshes edges
      auto matches = query::rtree::getClosestEdge(fVertices[i], searchSpace, nnearest);
      bool found   = false;
      for (const auto &match : matches) {
        auto weights = query::generateInterpolationElements(fVertices[i], tEdges[match.index]);
        if (std::all_of(weights.begin(), weights.end(), [](query::InterpolationElement const &elem) { return elem.weight >= 0.0; })) {
          _weights[i] = std::move(weights);
          distanceStatistics(match.distance);
          found = true;
          break;
        }
      }

      if (not found) {
        // Search for the origin inside the destination meshes vertices
        auto matchedVertices = query::rtree::getClosestVertex(fVertices[i], searchSpace);
        _weights[i]          = query::generateInterpolationElements(fVertices[i], tVertices[matchedVertices.at(0).index]);
        distanceStatistics(matchedVertices.at(0).distance);
      }
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  } else {
    const auto &tTriangles = searchSpace->triangles();
    if (!fVertices.empty() && tTriangles.empty()) {
      PRECICE_WARN("3D Mesh \"" << searchSpace->getName() << "\" does not contain triangles. Nearest projection mapping will map to primitives of lower dimension.");
    }

    // Lazy evaluation of indices for edges and vertices.
    // These are not necessary in the case of matching meshes.
    utils::statistics::DistanceAccumulator distanceStatistics;

    for (size_t i = 0; i < fVertices.size(); i++) {
      // Search for the vertex inside the destination meshes triangles
      auto matches = query::rtree::getClosestTriangle(fVertices[i], searchSpace, nnearest);
      bool found   = false;
      for (const auto &match : matches) {
        auto weights = query::generateInterpolationElements(fVertices[i], tTriangles[match.index]);
        if (std::all_of(weights.begin(), weights.end(), [](query::InterpolationElement const &elem) { return elem.weight >= 0.0; })) {
          _weights[i] = std::move(weights);
          found       = true;
          distanceStatistics(match.distance);
          break;
        }
      }

      if (not found) {
        // Search for the vertex inside the destination meshes edges
        matches.clear();
        matches = query::rtree::getClosestEdge(fVertices[i], searchSpace, nnearest);
        for (const auto &match : matches) {
          auto weights = query::generateInterpolationElements(fVertices[i], tEdges[match.index]);
          if (std::all_of(weights.begin(), weights.end(), [](query::InterpolationElement const &elem) { return elem.weight >= 0.0; })) {
            _weights[i] = std::move(weights);
            found       = true;
            distanceStatistics(match.distance);
            break;
          }
        }
      }

      if (not found) {
        // Search for the vertex inside the destination meshes vertices
        auto matchedVertices = query::rtree::getClosestVertex(fVertices[i], searchSpace);
        _weights[i]          = query::generateInterpolationElements(fVertices[i], tVertices[matchedVertices.at(0).index]);
        distanceStatistics(matchedVertices.at(0).distance);
      }
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  }
  _hasComputedMapping = true;
}

bool NearestProjectionMapping::hasComputedMapping() const
{
  return _hasComputedMapping;
}

void NearestProjectionMapping::clear()
{
  PRECICE_TRACE();
  _weights.clear();
  _hasComputedMapping = false;
}

void NearestProjectionMapping::map(
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

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Map consistent");
    PRECICE_ASSERT(_weights.size() == output()->vertices().size(),
                   _weights.size(), output()->vertices().size());
    for (size_t i = 0; i < output()->vertices().size(); i++) {
      InterpolationElements &elems     = _weights[i];
      size_t                 outOffset = i * dimensions;
      for (query::InterpolationElement &elem : elems) {
        size_t inOffset = (size_t) elem.element->getID() * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
          PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    PRECICE_DEBUG("Map conservative");
    PRECICE_ASSERT(_weights.size() == input()->vertices().size(),
                   _weights.size(), input()->vertices().size());
    for (size_t i = 0; i < input()->vertices().size(); i++) {
      size_t                 inOffset = i * dimensions;
      InterpolationElements &elems    = _weights[i];
      for (query::InterpolationElement &elem : elems) {
        size_t outOffset = (size_t) elem.element->getID() * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          PRECICE_ASSERT(outOffset + dim < (size_t) outValues.size());
          PRECICE_ASSERT(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  }
}

void NearestProjectionMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map.np.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  PRECICE_DEBUG("Compute Mapping for Tagging");

  computeMapping();
  PRECICE_DEBUG("Tagging First Round");

  // Determine the Mesh to Tag
  mesh::PtrMesh origins;
  if (getConstraint() == CONSISTENT) {
    origins = input();
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    origins = output();
  }

  // Gather all vertices to be tagged in a first phase.
  // max_count is used to shortcut if all vertices have been tagged.
  std::unordered_set<mesh::Vertex const *> tagged;
  const std::size_t                        max_count = origins->vertices().size();

  for (const InterpolationElements &elems : _weights) {
    for (const query::InterpolationElement &elem : elems) {
      if (!math::equals(elem.weight, 0.0)) {
        tagged.insert(elem.element);
      }
    }
    // Shortcut if all vertices are tagged
    if (tagged.size() == max_count) {
      break;
    }
  }

  // Now tag all vertices to be tagged in the second phase.
  for (auto &v : origins->vertices()) {
    if (tagged.count(&v) == 1) {
      v.tag();
    }
  }
  PRECICE_DEBUG("First Round Tagged " << tagged.size() << "/" << max_count << " Vertices");

  clear();
}

void NearestProjectionMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NP mapping no operation needed here
}

} // namespace mapping
} // namespace precice

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
#include "mesh/RTree.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "query/FindClosest.hpp"
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

namespace {
struct MatchType {
  double distance;
  int    index;

  MatchType() = default;
  MatchType(double d, int i)
      : distance(d), index(i){};
  constexpr bool operator<(MatchType const &other) const
  {
    return distance < other.distance;
  };
};
} // namespace

void NearestProjectionMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size(), output()->vertices().size());
  const std::string     baseEvent = "map.np.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  // Setup Direction of Mapping
  mesh::PtrMesh origins, search_space;
  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    origins      = output();
    search_space = input();
  } else {
    PRECICE_DEBUG("Compute conservative mapping");
    origins      = input();
    search_space = output();
  }

  const auto &fVertices = origins->vertices();
  const auto &tVertices = search_space->vertices();
  const auto &tEdges    = search_space->edges();

  _weights.resize(fVertices.size());

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  if (getDimensions() == 2) {
    if (!fVertices.empty() && tEdges.empty()) {
      PRECICE_WARN("2D Mesh \"" << search_space->getName() << "\" does not contain edges. Nearest projection mapping falls back to nearest neighbor mapping.");
    }

    precice::utils::Event e2(baseEvent + ".getIndexOnEdges", precice::syncMode);
    auto                  indexEdges = mesh::rtree::getEdgeRTree(search_space);
    e2.stop();

    // Lazy evaluation of the vertex index.
    // This is not necessary in the case of matching meshes.
    mesh::rtree::vertex_traits::Ptr indexVertices;

    utils::statistics::DistanceAccumulator distanceStatistics;

    std::vector<MatchType> matches;
    matches.reserve(nnearest);
    for (size_t i = 0; i < fVertices.size(); i++) {
      const Eigen::VectorXd &coords = fVertices[i].getCoords();
      // Search for the origin inside the destination meshes edges
      matches.clear();
      indexEdges->query(bg::index::nearest(coords, nnearest),
                        boost::make_function_output_iterator([&](int match) {
                          matches.emplace_back(bg::distance(coords, tEdges[match]), match);
                        }));
      std::sort(matches.begin(), matches.end());
      bool found = false;
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
        if (!indexVertices) {
          precice::utils::Event e3(baseEvent + ".getIndexOnVertices", precice::syncMode);
          indexVertices = mesh::rtree::getVertexRTree(search_space);
        }
        // Search for the origin inside the destination meshes vertices
        indexVertices->query(bg::index::nearest(coords, 1),
                             boost::make_function_output_iterator([&](int match) {
                               _weights[i] = query::generateInterpolationElements(fVertices[i], tVertices[match]);
                               distanceStatistics(bg::distance(fVertices[i], tVertices[match]));
                             }));
      }
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  } else {
    const auto &tTriangles = search_space->triangles();
    if (!fVertices.empty() && tTriangles.empty()) {
      PRECICE_WARN("3D Mesh \"" << search_space->getName() << "\" does not contain triangles. Nearest projection mapping will map to primitives of lower dimension.");
    }

    precice::utils::Event e2(baseEvent + ".getIndexOnTriangles", precice::syncMode);
    auto                  indexTriangles = mesh::rtree::getTriangleRTree(search_space);
    e2.stop();

    // Lazy evaluation of indices for edges and vertices.
    // These are not necessary in the case of matching meshes.
    mesh::rtree::edge_traits::Ptr   indexEdges;
    mesh::rtree::vertex_traits::Ptr indexVertices;

    utils::statistics::DistanceAccumulator distanceStatistics;

    std::vector<MatchType> matches;
    matches.reserve(nnearest);
    for (size_t i = 0; i < fVertices.size(); i++) {
      const Eigen::VectorXd &coords = fVertices[i].getCoords();

      // Search for the vertex inside the destination meshes triangles
      matches.clear();
      indexTriangles->query(bg::index::nearest(coords, nnearest),
                            boost::make_function_output_iterator([&](mesh::rtree::triangle_traits::IndexType const &match) {
                              matches.emplace_back(bg::distance(coords, tTriangles[match.second]), match.second);
                            }));
      std::sort(matches.begin(), matches.end());
      bool found = false;
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
        if (!indexEdges) {
          precice::utils::Event e3(baseEvent + ".getIndexOnEdges", precice::syncMode);
          indexEdges = mesh::rtree::getEdgeRTree(search_space);
        }
        // Search for the vertex inside the destination meshes edges
        matches.clear();
        indexEdges->query(bg::index::nearest(coords, nnearest),
                          boost::make_function_output_iterator([&](int match) {
                            matches.emplace_back(bg::distance(coords, tEdges[match]), match);
                          }));
        std::sort(matches.begin(), matches.end());
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
        if (!indexVertices) {
          precice::utils::Event e4(baseEvent + ".getIndexOnVertices", precice::syncMode);
          indexVertices = mesh::rtree::getVertexRTree(search_space);
        }
        // Search for the vertex inside the destination meshes vertices
        indexVertices->query(bg::index::nearest(coords, 1),
                             boost::make_function_output_iterator([&](int match) {
                               _weights[i] = query::generateInterpolationElements(fVertices[i], tVertices[match]);
                               distanceStatistics(bg::distance(fVertices[i], tVertices[match]));
                             }));
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

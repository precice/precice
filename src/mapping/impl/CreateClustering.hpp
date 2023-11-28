#pragma once

#include <Eigen/Core>

#include <algorithm>

#include "math/math.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "query/Index.hpp"
// required for the squared distance computation
#include "mapping/RadialBasisFctSolver.hpp"

namespace precice {
namespace mapping {
namespace impl {

using Vertices = std::vector<mesh::Vertex>;

/**
 * This file contains helper functions for the clustering of a given mesh as required for the
 * partition of unity data mapping class.
 */

// Internal helper functions
namespace {

/**
 * @brief Transforms a multi-dimensional index into a unique index. The formula here corresponds to a
 * dimension-major scheme. The index is used to put all generated cluster centers into a unique and
 * and accessible place. The following functions @ref zCurveNeighborOffsets also allow to immediately
 * access (spatially) adjacent clusters.
 *
 * @param[in] ids The multi-dimensional index of the cluster center (x, y, z)
 * @param[in] nCells The maximum size of each index (x_max, y_max, z_max)
 * @return VertexID The serialized unique ID
 */
constexpr VertexID zCurve(std::array<unsigned int, 3> ids, std::array<unsigned int, 3> nCells)
{
  return (ids[0] * nCells[1] + ids[1]) * nCells[2] + ids[2];
}

// Overload to handle the offsets (via ints) properly
constexpr VertexID zCurve(std::array<int, 3> ids, std::array<unsigned int, 3> nCells)
{
  return (ids[0] * nCells[1] + ids[1]) * nCells[2] + ids[2];
}

/**
 * @brief Computes the index offsets of all spatially adjacent cluster centers
 *
 * @tparam dim spatial dimension of the multi-dimensional index
 * @param nCells the maximum size of each index (x_max, y_max, z_max)
 * @return The serialized unique neighbor ids
 */
template <int dim>
std::array<int, math::pow_int<dim>(3) - 1> zCurveNeighborOffsets(std::array<unsigned int, 3> nCells);

/// 2D implementation
template <>
std::array<int, 8> zCurveNeighborOffsets<2>(std::array<unsigned int, 3> nCells)
{
  // All  neighbor directions
  const std::array<std::array<int, 3>, 8> neighbors2DIndex = {{{-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, -1, 0}, {0, 1, 0}, {1, -1, 0}, {1, 0, 0}, {1, 1, 0}}};
  std::array<int, 8>                      neighbors2D{{}};
  // Uses the int overload of the zCurve
  // @todo: mark the function as constexpr, when moving to C++20 (transform is non-constexpr)
  std::transform(neighbors2DIndex.begin(), neighbors2DIndex.end(), neighbors2D.begin(), [&nCells](auto &v) { return zCurve(v, nCells); });

  return neighbors2D;
}

/// 3D implementation
template <>
std::array<int, 26> zCurveNeighborOffsets<3>(std::array<unsigned int, 3> nCells)
{
  // All  neighbor directions
  const std::array<std::array<int, 3>, 26> neighbors3DIndex = {{{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 0, -1}, {-1, 1, 0}, {-1, 1, 1}, {-1, 1, -1}, {0, -1, -1}, {0, -1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, 1, 1}, {0, 1, -1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 0, -1}, {1, 1, 0}, {1, 1, 1}, {1, 1, -1}}};
  std::array<int, 26>                      neighbors3D{{}};
  // Uses the int overload of the zCurve
  std::transform(neighbors3DIndex.begin(), neighbors3DIndex.end(), neighbors3D.begin(), [&nCells](auto &v) { return zCurve(v, nCells); });

  return neighbors3D;
}
/**
 * @brief This function computes all vertices (in \p vertices) within a sphere. The sphere is defined through a center vertex ( \p centerID ),
 * which needs to be part of the search space/vertex container itself and a \p radius . The function checks all vertices given in \p neighborOffsets.
 *
 * @tparam ArrayType static array of different sizes, depending on the spatial dimension
 *
 * @param[in] vertices The vertex container containing the center vertex ( \p centerID ) and the search space
 * @param[in] centerID ID of the center vertex defining the sphere
 * @param[in] radius radius of the sphere
 * @param[in] neighborOffsets Index offsets of adjacent vertices (depends only on the spatial dimension, see also \ref zCurveNeighborOffsets())
 *
 * @return std::vector<VertexID> a collection of vertexIDs (referring again to \p vertices ) of vertices within the sphere.
 */
template <typename ArrayType>
std::vector<VertexID> getNeighborsWithinSphere(const Vertices &vertices, VertexID centerID, double radius, const ArrayType neighborOffsets)
{
  // number of neighbors = (3^dim)-1
  static_assert((neighborOffsets.size() == 26) || (neighborOffsets.size() == 8));
  std::vector<VertexID> result;
  for (auto off : neighborOffsets) {
    // corner case where we have 'dead axis' in the bounding box
    if (off == 0) {
      continue;
    }
    // compute neighbor ID
    // @todo target type should be an unsigned int
    PRECICE_ASSERT(off != 0);
    // Assumption, vertex position in the container == vertex ID
    VertexID neighborID = centerID + off;
    // We might run out of the array bounds along the edges, so we only consider vertices inside
    if (neighborID >= 0 && neighborID < static_cast<int>(vertices.size())) {
      auto &neighborVertex = vertices[neighborID];
      auto &center         = vertices[centerID];
      if (!neighborVertex.isTagged() && (computeSquaredDifference(center.rawCoords(), neighborVertex.rawCoords()) < math::pow_int<2>(radius))) {
        result.emplace_back(neighborID);
      }
    }
  }
  return result;
}

/**
 * @brief Given vertices in the \p clusterCenters and a \p clusterRadius , this function
 * tags vertices in the \p clusterCenters, which don't contain any vertex of the \p mesh
 *
 * @param[in] clusterCenters vertices representing the cluster centers
 * @param[in] clusterRadius radius of the clusters
 * @param[in] mesh Mesh, in which we look for vertices
 */
void tagEmptyClusters(Vertices &clusterCenters, double clusterRadius, mesh::PtrMesh mesh)
{
  // Alternative implementation: mesh->index().getVerticesInsideBox() == 0
  std::for_each(clusterCenters.begin(), clusterCenters.end(), [&](auto &v) {
    if (!v.isTagged() && !mesh->index().isAnyVertexInsideBox(v, clusterRadius)) {
      v.tag();
    }
  });
}

/**
 * @brief Projects non-tagged cluster centers to the closest vertices of the given \p mesh
 *
 * @param[in] clusterCenters the collection of cluster centers we want to project
 * @param[in] mesh the mesh we search for the closest vertex in order to move the center
 *
 * @note This function preserves the size of the \p clusterCenters , as tagged vertices remain unchanged.
 */
void projectClusterCentersToinputMesh(Vertices &clusterCenters, mesh::PtrMesh mesh)
{
  std::transform(clusterCenters.begin(), clusterCenters.end(), clusterCenters.begin(), [&](auto &v) {
    if (!v.isTagged()) {
      auto closestCenter = mesh->index().getClosestVertex(v.getCoords()).index;
      return mesh::Vertex{mesh->vertices()[closestCenter].getCoords(), v.getID()};
    } else {
      return v;
    }
  });
}

/**
 * @brief Tag (non-tagged) duplicated vertices in the \p centers , where duplicate is measured as a distance
 * between vertices smaller than a \p threshold . Here, we use the static neighbors from our index curve
 * (see also \ref zCurve and \ref zCurveNeighborOffsets ).
 *
 * @tparam dim spatial dimension of the clustering
 *
 * @param[in] centers vertex container, where we look for duplicate vertices
 * @param[in] nCells the maximum size of each index (x_max, y_max, z_max)
 * @param[in] threshold threshold value, which compares against the distance between centers
 */
template <int dim>
void tagDuplicateCenters(Vertices &centers, std::array<unsigned int, 3> nCells, double threshold)
{
  PRECICE_ASSERT(threshold >= 0);
  if (centers.empty())
    return;

  // Get the index offsets for the z curve
  auto neighborOffsets = zCurveNeighborOffsets<dim>(nCells);
  static_assert((neighborOffsets.size() == 8 && dim == 2) || (neighborOffsets.size() == 26 && dim == 3));

  //For the following to work, the vertexID has to correspond to the position in the centers
  PRECICE_ASSERT(std::all_of(centers.begin(), centers.end(), [idx = 0](auto &v) mutable { return v.getID() == idx++; }));
  // we check all neighbors
  for (auto &v : centers) {
    if (!v.isTagged()) {
      auto ids = getNeighborsWithinSphere(centers, v.getID(), threshold, neighborOffsets);
      // @todo check more selective which one to remove
      if (ids.size() > 0)
        v.tag();
    }
  }
}

/**
 * @brief Removes and erases tagged vertices from the input \p container.
 *
 * @param[in] container vertex container
 *
 * @note The neighbor search as carried out above (see \ref tagDuplicateCenters) won't work after calling this function,
 * as the position of vertices in the container changes.
 */
void removeTaggedVertices(Vertices &container)
{
  container.erase(std::remove_if(container.begin(), container.end(), [](auto &v) { return v.isTagged(); }), container.end());
}
} // namespace

/**
 * @brief Computes an estimate for the cluster radius, which results in approximately \p verticesPerCluster vertices inside
 * of each cluster. The algorithm generates random samples in the domain and queries the \p verticesPerCluster nearest-neighbors
 * from the mesh index tree. The cluster radius is then estimated through the distance between the center vertex (random sample)
 * and the vertex the furthest away from the center (being on the edge of the cluster).
 *
 * @param[in] verticesPerCluster target number of vertices in each cluster.
 * @param[in] inMesh mesh we want to create the clustering on
 * @param[in] bb bounding box of the domain. Used to place the random samples in the domain.
 *
 * @return The estimate for the cluster radius.
 */
inline double estimateClusterRadius(unsigned int verticesPerCluster, mesh::PtrMesh inMesh, const mesh::BoundingBox &bb)
{
  // Step 1: Generate random samples from the input mesh
  std::vector<VertexID> randomSamples;

  PRECICE_ASSERT(!bb.empty(), "Center computation won't work right now.");
  // All samples are 'moved' to the closest vertex from the input mesh in order
  // to avoid empty clusters when we have shell-like meshes (as in surface coupling)
  // The first sample: the vertex closest to the bounding box center
  const auto bbCenter = bb.center();
  randomSamples.emplace_back(inMesh->index().getClosestVertex(bbCenter).index);
  // Now we place samples for each space dimension above and below the bounding box center
  // and look for the closest vertex in the input mesh
  // In 2D the sampling would like this
  // +---------------------+
  // |          x          |
  // |                     |
  // |    x     c      x   |
  // |                     |
  // |          x          |
  // +---------------------+

  for (int d = 0; d < inMesh->getDimensions(); ++d) {
    auto sample = bbCenter;
    // Lower as the center
    sample[d] -= bb.getEdgeLength(d) * 0.25;
    randomSamples.emplace_back(inMesh->index().getClosestVertex(sample).index);
    // Higher than the center
    sample[d] += bb.getEdgeLength(d) * 0.5;
    randomSamples.emplace_back(inMesh->index().getClosestVertex(sample).index);
  }

  // Step 2: Compute the radius of the randomSamples ('centers'), which would have verticesPerCluster vertices
  std::vector<double> sampledClusterRadii;
  for (auto s : randomSamples) {
    // ask the index tree for the k-nearest neighbors  in order to estimate the point density
    auto kNearestVertexIDs = inMesh->index().getClosestVertices(inMesh->vertices()[s].getCoords(), verticesPerCluster);
    // compute the distance of each point to the center
    std::vector<double> squaredRadius(kNearestVertexIDs.size());
    std::transform(kNearestVertexIDs.begin(), kNearestVertexIDs.end(), squaredRadius.begin(), [&inMesh, s](auto i) {
      return computeSquaredDifference(inMesh->vertices()[i].rawCoords(), inMesh->vertices()[s].rawCoords());
    });
    // Store the maximum distance
    auto maxRadius = std::max_element(squaredRadius.begin(), squaredRadius.end());
    sampledClusterRadii.emplace_back(std::sqrt(*maxRadius));
  }

  // Step 3: Sort (using nth_element) the sampled radii and select the median as cluster radius
  PRECICE_ASSERT(sampledClusterRadii.size() % 2 != 0, "Median calculation is only valid for odd number of elements.");
  PRECICE_ASSERT(sampledClusterRadii.size() > 0);
  unsigned int middle = sampledClusterRadii.size() / 2;
  std::nth_element(sampledClusterRadii.begin(), sampledClusterRadii.begin() + middle, sampledClusterRadii.end());
  double clusterRadius = sampledClusterRadii[middle];

  return clusterRadius;
}

/**
 * @brief Creates a clustering as a collection of Vertices (representing the cluster centers) and a cluster radius,
 * as required for the partition of unity mapping. The algorithm estimates a cluster radius based on the input parameter
 * \p verticesPerCluster (see also \ref estimateClusterRadius above, which is directly used by the function). Afterwards,
 * the algorithm creates a cartesian-like grid of center vertices, where the distance of the centers is defined through
 * the \p relativeOverlap and the cluster radius. The parameter \p projectClustersToInput moves the cartesian center
 * vertices to the closest vertex from the input mesh, which is useful in case of very irregular meshes or shell-shaped
 * meshes.
 * The algorithm also removes potentially empty cluster, i.e., clusters which would have either no vertex from the
 * \p inMesh or from the \p outMesh . See also \ref tagEmptyClusters.
 *
 * @param[in] inMesh The input mesh (input mesh for consistent, output mesh for conservative mappings), on which the
 *            clustering is computed. The input parameters \p verticesPerCluster and \p projectClustersToInput refer
 *            to the \p inMesh
 * @param[in] outMesh The output mesh (output mesh for consistent, input mesh for conservative mappings),
 * @param[in] relativeOverlap Value between zero and one, which steers the relative distance between cluster centers.
 *            A value of zero leads to no overlap, a value of one would lead to a complete overlap between clusters.
 * @param[in] verticesPerCluster Target number of vertices per partition.
 * @param[in] projectClustersToInput if enabled, moves the cluster centers to the closest vertex of the \p inMesh
 *
 * @return a tuple for the cluster radius and a vector of vertices marking the cluster centers
 */
inline std::tuple<double, Vertices> createClustering(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh,
                                                     double relativeOverlap, unsigned int verticesPerCluster,
                                                     bool projectClustersToInput)
{
  precice::logging::Logger _log{"impl::createClustering"};
  PRECICE_TRACE();
  PRECICE_ASSERT(relativeOverlap < 1);
  PRECICE_ASSERT(verticesPerCluster > 0);
  PRECICE_ASSERT(inMesh->getDimensions() == outMesh->getDimensions());

  // If we have either no input or no output vertices, we return immediately
  if (inMesh->vertices().size() == 0 || outMesh->vertices().size() == 0)
    return {double{}, Vertices{}};

  PRECICE_ASSERT(!outMesh->vertices().empty() && !inMesh->vertices().empty());

  // startGridAtEdge boolean switch in order to decide either to start the clustering at the edge of the bounding box in each direction
  // (true) or start the clustering inside the bounding box (edge + 0.5 radius). The latter approach leads to fewer clusters,
  // but might result in a worse clustering
  const bool startGridAtEdge = projectClustersToInput;

  PRECICE_DEBUG("Relative overlap: {}", relativeOverlap);
  PRECICE_DEBUG("Vertices per cluster: {}", verticesPerCluster);

  // Step 1: Compute the local bounding box of the input mesh manually
  // Note that we don't use the corresponding bounding box functions from
  // precice::mesh (e.g. ::getBoundingBox), as the stored bounding box might
  // have the wrong size (e.g. direct access)
  // @todo: Which mesh should be used in order to determine the cluster centers:
  // pro outMesh: we want perfectly fitting clusters around our output vertices
  // however, this makes the cluster distribution/mapping dependent on the output
  precice::mesh::BoundingBox localBB = inMesh->index().getRtreeBounds();

#ifndef NDEBUG
  // Safety check
  precice::mesh::BoundingBox bb_check(inMesh->getDimensions());
  for (const mesh::Vertex &vertex : inMesh->vertices()) {
    bb_check.expandBy(vertex);
  }
  PRECICE_ASSERT(bb_check == localBB);
#endif

  // If we have very few vertices in the domain, (in this case twice our cluster
  // size as we decompose most probably at least in 4 clusters) we just use a
  // single cluster. The clustering result of the algorithm further down is in
  // this case not optimal and might lead to too many clusters.
  // The single cluster has in principle a radius of inf. We use here twice the
  // length of the longest bounding box edge length and the center of the bounding
  // box for the center point.
  if (inMesh->vertices().size() < verticesPerCluster * 2)
    return {localBB.longestEdgeLength() * 2, Vertices{mesh::Vertex({localBB.center(), 0})}};

  // We define a convenience alias for the localBB. In case we need to synchronize the clustering across ranks later on, we need
  // to work with the global bounding box of the whole domain.
  auto globalBB = localBB;

  // Step 2: Now we pick random samples from the input mesh and ask the index tree for the k-nearest neighbors
  // in order to estimate the point density and determine a proper cluster radius (see also the function documentation)
  double clusterRadius = estimateClusterRadius(verticesPerCluster, inMesh, localBB);
  PRECICE_DEBUG("Vertex cluster radius: {}", clusterRadius);

  // maximum distance between cluster centers lying diagonal to each other. The maximum distance takes the overlap condition into
  // account: if the distance between the centers is sqrt(2) * radius, we violate the overlap condition between diagonal clusters
  // 0.3 should be a good default value
  // @todo: most of the 3D experiments ran with the default of 0.3 which would correspond to relativeOverlap of 0, i.e., the limit
  // Can we use different values for 2D and 3D? Which value should we use?
  const double maximumCenterDistance = inMesh->getDimensions() == 2 ? std::sqrt(2) * clusterRadius * (1 - relativeOverlap) : clusterRadius * (1 - relativeOverlap);

  // Step 3: using the maximum distance and the bounding box, compute the number of clusters in each direction
  // we ceil the number of clusters in order to guarantee the desired overlap
  std::array<unsigned int, 3> nClustersGlobal{1, 1, 1};
  for (int d = 0; d < globalBB.getDimension(); ++d)
    nClustersGlobal[d] = std::ceil(std::max(1., globalBB.getEdgeLength(d) / maximumCenterDistance));

  // Step 4: Determine the centers of the clusters
  Vertices centers;
  // Vector used to temporarily store each center coordinates
  std::vector<double> centerCoords(inMesh->getDimensions());
  // Vector storing the distances between the cluster in each direction
  std::vector<double> distances(inMesh->getDimensions());
  // Vector storing the starting coordinates in each direction
  std::vector<double> start(inMesh->getDimensions());
  // Fill the constant vectos, i.e., the distances and the start
  for (unsigned int d = 0; d < distances.size(); ++d) {
    // this distance calculation guarantees that we have at least the minimal specified overlap, as we ceil the division for the number of clusters
    // as an alternative, once could use distance = const = maximumCenterDistance, which might lead to better results when projecting and removing duplicates (?)
    distances[d] = globalBB.getEdgeLength(d) / (nClustersGlobal[d]);
    start[d]     = globalBB.minCorner()(d);
  }
  PRECICE_DEBUG("Distances: {}", distances);

  // Step 5: Take care of the starting layer: if we start the grid at the globalBB edge we have an additional layer in each direction
  if (startGridAtEdge) {
    for (int d = 0; d < inMesh->getDimensions(); ++d) {
      if (globalBB.getEdgeLength(d) > math::NUMERICAL_ZERO_DIFFERENCE) {
        nClustersGlobal[d] += 1;
      }
    }
  } // If we don't start at the edge, we shift the starting layer into the middle by 0.5 * distance
  else {
    std::transform(start.begin(), start.end(), distances.begin(), start.begin(), [](auto s, auto d) { return s + 0.5 * d; });
  }

  // Print some information
  PRECICE_DEBUG("Global cluster distribution: {}", nClustersGlobal);
  unsigned int nTotalClustersGlobal = std::accumulate(nClustersGlobal.begin(), nClustersGlobal.end(), 1U, std::multiplies<unsigned int>());
  PRECICE_DEBUG("Global number of total clusters (tentative): {}", nTotalClustersGlobal);

  // Step 6: transform the global metrics (number of clusters and starting layer) into local ones
  // Since the local and global bounding boxes are the same in the current implementation, the
  // global metrics and the local metrics will be the same.
  // First, we need to update the starting points of our clustering scheme
  std::array<unsigned int, 3> nClustersLocal{1, 1, 1};
  for (unsigned int d = 0; d < start.size(); ++d) {
    // Exclude the case where we have no distances (nClustersGlobal[d] == 1 )
    if (distances[d] > 0) {
      start[d] += std::ceil(((localBB.minCorner()(d) - start[d]) / distances[d]) - math::NUMERICAL_ZERO_DIFFERENCE) * distances[d];
      // One cluster for the starting point and a further one for each distance we can fit into the BB
      nClustersLocal[d] = 1 + std::floor(((localBB.maxCorner()(d) - start[d]) / distances[d]) + math::NUMERICAL_ZERO_DIFFERENCE);
    }
  }

  unsigned int nTotalClustersLocal = std::accumulate(nClustersLocal.begin(), nClustersLocal.end(), 1U, std::multiplies<unsigned int>());
  centers.resize(nTotalClustersLocal, mesh::Vertex({centerCoords, -1}));

  // Print some information
  PRECICE_DEBUG("Local cluster distribution: {}", nClustersLocal);
  PRECICE_DEBUG("Local number of total clusters (tentative): {}", nTotalClustersLocal);

  // Step 7: fill the center container using the zCurve for indexing. For the further processing, it is also important to ensure
  // that the ID within the center Vertex class coincides with the position of the center vertex in the container.
  // We start with the (bottom left) corner and iterate over all dimension
  if (inMesh->getDimensions() == 2) {
    centerCoords[0] = start[0];
    for (unsigned int x = 0; x < nClustersLocal[0]; ++x, centerCoords[0] += distances[0]) {
      centerCoords[1] = start[1];
      for (unsigned int y = 0; y < nClustersLocal[1]; ++y, centerCoords[1] += distances[1]) {
        auto id = zCurve(std::array<unsigned int, 3>{{x, y, 0}}, nClustersLocal);
        PRECICE_ASSERT(id < static_cast<int>(centers.size()) && id >= 0);
        centers[id] = mesh::Vertex({centerCoords, id});
      }
    }
  } else {
    centerCoords[0] = start[0];
    for (unsigned int x = 0; x < nClustersLocal[0]; ++x, centerCoords[0] += distances[0]) {
      centerCoords[1] = start[1];
      for (unsigned int y = 0; y < nClustersLocal[1]; ++y, centerCoords[1] += distances[1]) {
        centerCoords[2] = start[2];
        for (unsigned int z = 0; z < nClustersLocal[2]; ++z, centerCoords[2] += distances[2]) {
          auto id = zCurve(std::array<unsigned int, 3>{{x, y, z}}, nClustersLocal);
          PRECICE_ASSERT(id < static_cast<int>(centers.size()) && id >= 0);
          centers[id] = mesh::Vertex({centerCoords, id});
        }
      }
    }
  }

  // Step 8: tag vertices, which should be removed later on

  // Step 8a: tag empty clusters
  tagEmptyClusters(centers, clusterRadius, inMesh);
  tagEmptyClusters(centers, clusterRadius, outMesh);
  PRECICE_DEBUG("Number of non-tagged centers after empty clusters were filtered: {}", std::count_if(centers.begin(), centers.end(), [](auto &v) { return !v.isTagged(); }));

  // Step 9: Move the cluster centers if desired
  if (projectClustersToInput) {
    projectClusterCentersToinputMesh(centers, inMesh);
    // @todo: The duplication should probably be defined in terms of the target distance, not the actual distance. Find good default values
    const auto duplicateThreshold = 0.4 * (*std::min_element(distances.begin(), distances.end()));
    PRECICE_DEBUG("Tagging duplicates using the threshold value {} for the regular distances {}", duplicateThreshold, distances);
    // usually inMesh->getDimensions() == 2
    // to also cover the case where we have only a single layer in a 3D computation we use the number of clusters in z direction
    // Step 8b or 9a: Moving cluster centers might lead to duplicate centers or centers being too close to each other. Therefore,
    // we tag duplicate centers (see also the documentation of \ref tagDuplicateCenters ).
    if (nClustersLocal[2] == 1) {
      tagDuplicateCenters<2>(centers, nClustersLocal, duplicateThreshold);
    } else {
      tagDuplicateCenters<3>(centers, nClustersLocal, duplicateThreshold);
    }
    // the tagging here won't filter out a lot of clusters, but finding them anyway allows us to allocate the right amount of
    // memory for the cluster vector in the mapping, which would otherwise be expensive
    // we cannot hit any empty vertices in the inputmesh
    tagEmptyClusters(centers, clusterRadius, outMesh);
    PRECICE_DEBUG("Number of non-tagged centers after duplicate tagging : {}", std::count_if(centers.begin(), centers.end(), [](auto &v) { return !v.isTagged(); }));
  }
  PRECICE_ASSERT(std::all_of(centers.begin(), centers.end(), [idx = 0](auto &v) mutable { return v.getID() == idx++; }));

  // Step 10: finally, remove the vertices from the center container
  removeTaggedVertices(centers);
  PRECICE_ASSERT(std::none_of(centers.begin(), centers.end(), [](auto &v) { return v.isTagged(); }));
  PRECICE_CHECK(centers.size() > 0, "Too many vertices have been filtered out.");

  return {clusterRadius, centers};
}
} // namespace impl
} // namespace mapping
} // namespace precice

#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include "math/math.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "query/Index.hpp"
// required for the squared distance computation
#include "mapping/RadialBasisFctSolver.hpp"

namespace precice {
namespace mapping {
namespace partitioner {

using Vertices = std::vector<mesh::Vertex>;
namespace {
precice::logging::Logger _log{"partitioner::createUniformBlockPartitioning"};

// Formula doesn't correspond to a z curve right now, but it works
constexpr VertexID zCurve(std::array<unsigned int, 3> ids, std::array<unsigned int, 3> nCells)
{
  return (ids[0] * nCells[1] + ids[1]) * nCells[2] + ids[2];
}
// Overload to handle the offsets (via ints) properly
constexpr int zCurve(std::array<int, 3> ids, std::array<unsigned int, 3> nCells)
{
  return (ids[0] * nCells[1] + ids[1]) * nCells[2] + ids[2];
}

template <int dim>
std::array<int, math::pow_int<dim>(3) - 1> zCurveNeighborOffsets(std::array<unsigned int, 3> nCells);

template <>
std::array<int, 8> zCurveNeighborOffsets<2>(std::array<unsigned int, 3> nCells)
{
  // All  neighbor directions
  const std::array<std::array<int, 3>, 8> neighbors2DIndex = {{{-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, -1, 0}, {0, 1, 0}, {1, -1, 0}, {1, 0, 0}, {1, 1, 0}}};
  std::array<int, 8>                      neighbors2D{{}};
  // Uses the int overload of the zCurve
  std::transform(neighbors2DIndex.begin(), neighbors2DIndex.end(), neighbors2D.begin(), [&nCells](auto &v) { return zCurve(v, nCells); });

  return neighbors2D;
}

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

// Assumption, vertex position in the container == vertex ID
template <typename ArrayType>
std::vector<VertexID> getNeighborsWithinSphere(const Vertices &vertices, mesh::Vertex &center, const ArrayType &neighborOffsets, double radius)
{
  // number of neighbors = (3^dim)-1
  PRECICE_ASSERT((neighborOffsets.size() == 26) || (neighborOffsets.size() == 8));
  std::vector<VertexID> result;
  for (auto off : neighborOffsets) {
    // corner case where we have 'dead axis' in the bounding box
    if (off == 0) {
      continue;
    }
    // compute neighbor ID
    // @todo target type should be an unsigned int
    PRECICE_ASSERT(off != 0);
    VertexID neighborID = center.getID() + off;
    PRECICE_ASSERT(center.getID() == vertices[center.getID()].getID());
    // We might run out of the array bounds along the edges, so we only consider vertices inside
    if (neighborID >= 0 && neighborID < vertices.size()) {
      auto &neighborVertex = vertices[neighborID];
      if (!neighborVertex.isTagged() && (computeSquaredDifference(center.rawCoords(), neighborVertex.rawCoords()) < math::pow_int<2>(radius))) {
        result.emplace_back(neighborID);
      }
    }
  }
  return result;
}

/**
 * @brief Tags empty partitions, which are defined as no output vertices and 0.1size<targetSize
 */
void tagEmptyPartitions(double partitionRadius, Vertices &container, mesh::PtrMesh mesh, unsigned int threshold)
{
  PRECICE_TRACE();
  // todo: we don't need all vertices for the query, more than one would already allow to abort
  // todo: investigate difference to getClosestVertex > radius
  std::for_each(container.begin(), container.end(), [&](auto &v) {
    if (mesh->index().getVerticesInsideBox(v, partitionRadius).size() <= threshold) {
      v.tag();
    }
  });
}

/**
 * @brief Projects non-tagged partition centers to the closest vertices of the given mesh
 */
void projectPartitionCentersToinputMesh(Vertices &container, mesh::PtrMesh mesh)
{
  PRECICE_TRACE();
  std::transform(container.begin(), container.end(), container.begin(), [&](auto &v) {
    if (!v.isTagged()) {
      auto closestCenter = mesh->index().getClosestVertex(v.getCoords()).index;
      return mesh::Vertex{mesh->vertices()[closestCenter].getCoords(), v.getID()};
    } else {
      return v;
    }
  });
}

/**
 * @brief Tag non-tagged duplicated center vertices, where duplicate is measured as a distance smaller than a threshold
 * Here, we use the static neighbors from our partitioning curve
 */
template <int dim>
void tagDuplicateCenters(Vertices &container, double threshold, std::array<unsigned int, 3> nCells)
{
  PRECICE_TRACE();
  if (container.empty())
    return;

  // Get the index offsets for the z curve
  auto neighborOffsets = zCurveNeighborOffsets<dim>(nCells);
  PRECICE_ASSERT((neighborOffsets.size() == 8 && dim == 2) || (neighborOffsets.size() == 26 && dim == 3));

  // we check all neighbors
  for (auto &v : container) {
    if (!v.isTagged()) {
      auto ids = getNeighborsWithinSphere(container, v, neighborOffsets, threshold);
      // @todo check more selective which one to remove
      if (ids.size() > 0)
        v.tag();
    }
  }
}

/**
 * @brief Removes tagged vertices from the input container.
 * The neighbor search as carried out above won't work after calling this  function.
 *
 * @param container vertex container
 */
void removeTaggedVertices(Vertices &container)
{
  container.erase(std::remove_if(container.begin(), container.end(), [](auto &v) { return v.isTagged(); }), container.end());
}

/**
 * @brief Computes (via a bounding box broadcast) the global bounding box of the distributed mesh. The mesh
 *
 * @param mesh
 * @return mesh::BoundingBox
 */
mesh::BoundingBox getGlobalBoundingBox(mesh::PtrMesh mesh)
{
  PRECICE_ASSERT(!mesh->getBoundingBox().empty());

  // @TODO: We need a scatter from and to the primary rank of bounding boxes here

  return mesh->getBoundingBox();
}

double estimatePartitionRadius(unsigned int verticesPerPartition, mesh::PtrMesh inMesh, const mesh::BoundingBox &bb)
{
  // 2. Now we pick random samples from the input mesh and ask the index tree for the k-nearest neighbors
  // in order to estimate the point density and determine a proper partition radius
  std::vector<VertexID> randomSamples;

  PRECICE_ASSERT(!bb.empty(), "Center computation won't work right now.");
  // Sample the vertex closest to the bb center
  const auto bbCenter = bb.center();
  randomSamples.emplace_back(inMesh->index().getClosestVertex(bbCenter).index);
  for (int d = 0; d < inMesh->getDimensions(); ++d) {
    auto sample = bbCenter;
    // Lower as the center
    sample[d] -= bb.getEdgeLength(d) * 0.25;
    randomSamples.emplace_back(inMesh->index().getClosestVertex(sample).index);
    // Higher than the center
    sample[d] += bb.getEdgeLength(d) * 0.5;
    randomSamples.emplace_back(inMesh->index().getClosestVertex(sample).index);
  }

  // TODO: Add more samples here in order to get a better estimate

  // Compute the radius of the randomSample ('center'), which would have verticesPerPartition vertices
  std::vector<double> sampledPartitionRadii;
  for (auto s : randomSamples) {
    // Get the k-nearest neighbors
    auto kNearestVertexIDs = inMesh->index().getClosestVertices(inMesh->vertices()[s], verticesPerPartition);
    // Compute the distance of each point to the center
    std::vector<double> squaredRadius(kNearestVertexIDs.size());
    std::transform(kNearestVertexIDs.begin(), kNearestVertexIDs.end(), squaredRadius.begin(), [&inMesh, s](auto i) {
      return computeSquaredDifference(inMesh->vertices()[i].rawCoords(), inMesh->vertices()[s].rawCoords(), {{true, true, true}});
    });
    // Store the maximum distance
    auto maxRadius = std::max_element(squaredRadius.begin(), squaredRadius.end());
    sampledPartitionRadii.emplace_back(std::sqrt(*maxRadius));
  }

  // Select the median as partition radius
  PRECICE_ASSERT(sampledPartitionRadii.size() % 2 != 0, "Median calculation is only valid for odd number of elements.");
  PRECICE_ASSERT(sampledPartitionRadii.size() > 0);
  unsigned int middle = sampledPartitionRadii.size() / 2;
  std::nth_element(sampledPartitionRadii.begin(), sampledPartitionRadii.begin() + middle, sampledPartitionRadii.end());
  double averagePartitionRadius = sampledPartitionRadii[middle];

  PRECICE_INFO("Partition Radius: {}", averagePartitionRadius);
  return averagePartitionRadius;
}
} // namespace

/**
 * @brief Create a Uniform Block Partitioning using random smaples for the points density, given the meshes
 *
 * @param inMesh The input mesh (remote mesh for consistent and conservative)
 * @param outMesh The output mesh (local mesh for consistent and conservative)
 * @param projectPartitionsToInput if enabled, places the partition centers at the closest vertex of the input mesh
 *
 * @return a tuple for the partition radius and a vector of points marking the partition centers Vertices
 */
inline std::tuple<double, Vertices> createUniformBlockPartitioning(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh,
                                                                   double relativeOverlap, unsigned int verticesPerPartition,
                                                                   bool projectPartitionsToInput)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(relativeOverlap < 1);
  PRECICE_ASSERT(verticesPerPartition > 0);
  PRECICE_ASSERT(inMesh->getDimensions() == outMesh->getDimensions());

  if (inMesh->vertices().size() == 0 || outMesh->vertices().size() == 0)
    return {double{}, Vertices{}};

  PRECICE_ASSERT(!outMesh->vertices().empty() && !inMesh->vertices().empty());

  // startGridAtEdge boolean switch in order to decide either to start the partitioning at the edge of the bounding box in each direction
  // (true) or start the partitioning inside the bounding box (edge + 0.5 radius). The latter approach leads to fewer partitions,
  // but might result in a worse partitioning
  const bool startGridAtEdge = projectPartitionsToInput;

  // Get the number of points in the input mesh lying within the output mesh region
  // in order to estimate the global number of vertices this rank has to compute the
  // mapping on
  PRECICE_DEBUG("Relative overlap: {}", relativeOverlap);
  PRECICE_DEBUG("Vertices per partition: {}", verticesPerPartition);
  // 1. Get the (global) bounding box of the output mesh
  inMesh->computeBoundingBox();
  auto localBB = inMesh->getBoundingBox();

  // @todo: Which mesh should be used in order to determine the partition centers:
  // pro outMesh: we want perfectly fitting partitions around our output vertices
  // however, this makes the partition distribution/mapping dependent on the output
  auto globalBB = getGlobalBoundingBox(inMesh);

  // 2. Now we pick random samples from the input mesh and ask the index tree for the k-nearest neighbors
  // in order to estimate the point density and determine a proper partition radius
  double averagePartitionRadius = estimatePartitionRadius(verticesPerPartition, inMesh, localBB);

  // maximum distance between partition centers, which corresponds to the overlap condition, if the distance between the centers is sqrt(2) * radius,
  // we violate the overlap condition between diagonal partitions
  // 0.3 should be a good default value
  const double maximumCenterDistance = boost::math::constants::root_two<double>() * averagePartitionRadius * (1 - relativeOverlap);

  std::array<unsigned int, 3> nPartitionsGlobal{1, 1, 1};
  for (int d = 0; d < globalBB.getDimension(); ++d)
    nPartitionsGlobal[d] = std::ceil(std::max(1., globalBB.getEdgeLength(d) / maximumCenterDistance));

  // 1. Determine the centers of the partitions
  Vertices centers;
  // Store the center coordinates
  std::vector<double> centerCoords(inMesh->getDimensions());
  // Store the distance between each partition in each direction
  std::vector<double> distances(inMesh->getDimensions());
  // Store the starting point
  std::vector<double> start(inMesh->getDimensions());
  for (unsigned int d = 0; d < distances.size(); ++d) {
    // this distance calculation guarantees that we have at least the minimal specified overlap, as we ceil the division for the number of partitions
    // as an alternative, once could use distance = const = maximumCenterDistance, which might lead to better results when projecting and removing duplicates (?)
    distances[d] = globalBB.getEdgeLength(d) / (nPartitionsGlobal[d]);
    start[d]     = globalBB.getDirectionsCoordinates(d).first;
  }
  PRECICE_DEBUG("Distances: {}", distances);

  // If we start the grid at the globalBB edge we have an additional layer in each direction
  if (startGridAtEdge) {
    for (unsigned int d = 0; d < inMesh->getDimensions(); ++d) {
      if (globalBB.getEdgeLength(d) > math::NUMERICAL_ZERO_DIFFERENCE) {
        nPartitionsGlobal[d] += 1;
      }
    }
  } // If we don't start at the edge, we shift the starting layer into the middle by 0.5 * distance
  else {
    std::transform(start.begin(), start.end(), distances.begin(), start.begin(), [](auto s, auto d) { return s + 0.5 * d; });
  }

  // Print some information
  PRECICE_DEBUG("Global partition distribution: {}", nPartitionsGlobal);
  unsigned int nTotalPartitionsGlobal = std::accumulate(nPartitionsGlobal.begin(), nPartitionsGlobal.end(), 1U, std::multiplies<unsigned int>());
  PRECICE_DEBUG("Global number of total partitions (tentative): {}", nTotalPartitionsGlobal);

  // Now we need to transform the global metric into local ones
  // First, we need to update the starting points of our partitioning scheme
  PRECICE_DEBUG("Start before: {}", start);
  PRECICE_DEBUG("Distances: {}", distances);
  PRECICE_DEBUG("At edge? : {}", startGridAtEdge);
  std::array<unsigned int, 3> nPartitionsLocal{1, 1, 1};
  for (unsigned int d = 0; d < start.size(); ++d) {
    // Exclude the case where we have no distances (nPartitionsGlobal[d] == 1 )
    if (distances[d] > 0) {
      start[d] += std::ceil(((localBB.getDirectionsCoordinates(d).first - start[d]) / distances[d]) - math::NUMERICAL_ZERO_DIFFERENCE) * distances[d];
      // One partition for the starting point and a further one for each distance we can fit into the BB
      nPartitionsLocal[d] = 1 + std::floor(((localBB.getDirectionsCoordinates(d).second - start[d]) / distances[d]) + math::NUMERICAL_ZERO_DIFFERENCE);
    }
  }

  // Print some information
  PRECICE_DEBUG("Local partition distribution: {}", nPartitionsLocal);
  unsigned int nTotalPartitionsLocal = std::accumulate(nPartitionsLocal.begin(), nPartitionsLocal.end(), 1U, std::multiplies<unsigned int>());
  PRECICE_DEBUG("Local number of total partitions (tentative): {}", nTotalPartitionsLocal);

  centers.resize(nTotalPartitionsLocal, mesh::Vertex({centerCoords, -1}));

  // start with the (bottom left) corner
  if (inMesh->getDimensions() == 2) {
    centerCoords[0] = start[0];
    for (unsigned int x = 0; x < nPartitionsLocal[0]; ++x, centerCoords[0] += distances[0]) {
      centerCoords[1] = start[1];
      for (unsigned int y = 0; y < nPartitionsLocal[1]; ++y, centerCoords[1] += distances[1]) {
        auto id = zCurve(std::array<unsigned int, 3>{{x, y, 0}}, nPartitionsLocal);
        PRECICE_ASSERT(id < centers.size() && id >= 0);
        centers[id] = mesh::Vertex({centerCoords, id});
      }
    }
  } else {
    centerCoords[0] = start[0];
    for (unsigned int x = 0; x < nPartitionsLocal[0]; ++x, centerCoords[0] += distances[0]) {
      centerCoords[1] = start[1];
      for (unsigned int y = 0; y < nPartitionsLocal[1]; ++y, centerCoords[1] += distances[1]) {
        centerCoords[2] = start[2];
        for (unsigned int z = 0; z < nPartitionsLocal[2]; ++z, centerCoords[2] += distances[2]) {
          auto id = zCurve(std::array<unsigned int, 3>{{x, y, z}}, nPartitionsLocal);
          PRECICE_ASSERT(id < centers.size() && id >= 0);
          centers[id] = mesh::Vertex({centerCoords, id});
        }
      }
    }
  }

  tagEmptyPartitions(averagePartitionRadius, centers, inMesh, 0U);
  tagEmptyPartitions(averagePartitionRadius, centers, outMesh, 0U);
  PRECICE_DEBUG("Number of non-tagged centers after empty partitions were filtered: {}", std::count_if(centers.begin(), centers.end(), [](auto &v) { return !v.isTagged(); }));
  if (projectPartitionsToInput) {
    projectPartitionCentersToinputMesh(centers, inMesh);
    // @todo: The duplication should probably be defined in terms of the target distance, not the actual distance. Find good default values
    const auto duplicateThreshold = 0.4 * (*std::min_element(distances.begin(), distances.end()));
    PRECICE_DEBUG("Tagging duplicates using the threshold value {} fo the distances {}", duplicateThreshold, distances);
    // usually inMesh->getDimensions() == 2
    // to also cover the case where we have only a single layer in a 3D computation we use the number of partitions in z direction
    if (nPartitionsLocal[2] == 1) {
      tagDuplicateCenters<2>(centers, duplicateThreshold, nPartitionsLocal);
    } else {
      tagDuplicateCenters<3>(centers, duplicateThreshold, nPartitionsLocal);
    }
    PRECICE_DEBUG("Number of non-tagged centers after duplicate tagging : {}", std::count_if(centers.begin(), centers.end(), [](auto &v) { return !v.isTagged(); }));
  }
  PRECICE_ASSERT(std::all_of(centers.begin(), centers.end(), [idx = 0](auto &v) mutable { return v.getID() == idx++; }));
  removeTaggedVertices(centers);
  PRECICE_ASSERT(std::none_of(centers.begin(), centers.end(), [](auto &v) { return v.isTagged(); }));
  PRECICE_CHECK(centers.size() > 0, "Too many vertices have been filtered out.");

  return {averagePartitionRadius, centers};
}

/**
 * @brief Create a Uniform Block Partitioning in 2D, given the meshes
 *
 * @param inMesh The input mesh (remote mesh for consistent and conservative)
 * @param outMesh The output mesh (local mesh for consistent and conservative)
 *
 * @return a tuple for the partition radius and a vector of points marking the partition centers Vertices
 */
inline std::tuple<double, Vertices> createUniformBlockPartitioning2D(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh, double relativeOverlap, unsigned int verticesPerPartition)
{
  PRECICE_ASSERT(relativeOverlap < 1);
  PRECICE_ASSERT(verticesPerPartition > 0);
  PRECICE_ASSERT(inMesh->getDimensions() == outMesh->getDimensions());
  PRECICE_DEBUG("Creating uniform block partitioning");

  // Get the number of points in the input mesh lying within the output mesh region
  // in order to estimate the global number of vertices this rank has to compute the
  // mapping on
  PRECICE_DEBUG("Relative overlap: {}", relativeOverlap);
  PRECICE_DEBUG("Vertices per partition: {}", verticesPerPartition);
  // 1. Get the global bounding box of the output mesh
  outMesh->computeBoundingBox();
  auto               bb         = outMesh->getBoundingBox();
  auto               vertices   = inMesh->index().getVerticesInsideBox(bb);
  const unsigned int vertexSize = vertices.size();

  // 2. Based on the input parameter verticesPerPartition and overlap, estimate the number
  // of partitions and the radius
  PRECICE_INFO("Input mesh size: {}", vertexSize);
  unsigned int nTotalPartitions = std::max(1., (vertexSize / verticesPerPartition) * (1. / (1 - relativeOverlap)));
  PRECICE_INFO("Number of total partitions: {}", nTotalPartitions);

  // Get the edge length of the bounding box in each direction
  std::vector<double> edgeLength;
  for (unsigned int d = 0; d < bb.getDimension(); ++d) {
    edgeLength.emplace_back(bb.getEdgeLength(d));
  }
  PRECICE_ASSERT(bb.getDimension() == 2, "Not implemented.");

  // Assume uniform distribution in the bounding box:
  // xEdge/yEdge = nPartitionsX/nPartitionsY, nPartitionsX * nPartitionsY = nPartitions
  // --> nPartitionsX = sqrt(nPartitions* (xEdge/yEdge))
  std::array<unsigned int, 3> nPartitions{1, 1, 1};
  // Values are ceiled
  nPartitions[0] = std::ceil(std::max(1., std::sqrt(nTotalPartitions * (edgeLength[0] / edgeLength[1]))));
  nPartitions[1] = std::ceil(std::max(1., std::sqrt(nTotalPartitions * (edgeLength[1] / edgeLength[0]))));

  PRECICE_INFO("Partition distribution: {}", nPartitions);

  // Compute the radius based on the edge length and the number of partitions
  // 0.5 since we use a radius here instead of a diameter
  PRECICE_ASSERT(!bb.empty());
  const double distance               = std::max(edgeLength[0] / nPartitions[0], edgeLength[1] / nPartitions[1]);
  double       averagePartitionRadius = std::pow(relativeOverlap, 2) * distance + relativeOverlap * 0.5 * distance + 0.5 * distance;
  PRECICE_INFO("Partition Radius: {}", averagePartitionRadius);

  Vertices centers;
  {
    // 1. Determine the centers of the partitions
    const auto &        bb = outMesh->getBoundingBox();
    std::vector<double> centerCoords(inMesh->getDimensions());
    const double        center_x = bb.getEdgeLength(0) / nPartitions[0];
    const double        center_y = bb.getEdgeLength(1) / nPartitions[1];

    // start with the (bottom left) corner
    centerCoords[0] = bb.getDirectionsCoordinates(0).first + 0.5 * center_x;
    for (unsigned int x = 0; x < nPartitions[0]; ++x) {
      centerCoords[1] = bb.getDirectionsCoordinates(1).first + 0.5 * center_y;
      for (unsigned int y = 0; y < nPartitions[1]; ++y) {
        centers.emplace_back(mesh::Vertex({centerCoords, -1}));

        // for all x coordinates, iterate over the corresponding y coordinates
        // 2 since we are dealing with the radius from both partitions
        centerCoords[1] += center_y;
      }
      centerCoords[0] += center_x;
    }
  }
  return {averagePartitionRadius, centers};
}

} // namespace partitioner
} // namespace mapping
} // namespace precice

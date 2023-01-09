#pragma once

#include <Eigen/Core>
#include <numeric>

#include "com/Communication.hpp"
#include "io/ExportVTU.hpp"
#include "mapping/impl/CreateClustering.hpp"
#include "mapping/impl/SphericalVertexCluster.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/IntraComm.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

/**
 * Mapping using partition of unity decomposition strategies: The class here inherits from the Mapping
 * class and orchestrates the partitions (called vertex clusters) in order to represent a partition of unity.
 * This means in particular that the class computes the weights for the evaluation vertices and the necessary
 * association between evaluation vertices and the clustering.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class PartitionOfUnityMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimension Dimensionality of the meshes
   * @param[in] parameter shape parameter or support radius of the RB function
   * @param[in] verticesPerCluster Target number of vertices to be clustered together
   * @param[in] relativeOverlap Overlap between clusters, where 1 would correspond to a complete overlap, 0 to distance of 2 x radius between clusters
   */
  PartitionOfUnityMapping(
      Mapping::Constraint constraint,
      int                 dimension,
      double              parameter,
      Polynomial          polynomial,
      unsigned int        verticesPerCluster,
      double              relativeOverlap,
      bool                projectToInput);

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// tag the vertices required for the mapping
  virtual void tagMeshFirstRound() override;

  /// nothing to do here
  virtual void tagMeshSecondRound() override;

private:
  /// logger, as usual
  precice::logging::Logger _log{"mapping::PartitionOfUnityMapping"};

  /// main data container storing all the clusters, which need to be solved individually
  std::vector<SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>> _clusters;

  // Shape parameter or support radius for the RBF interpolant,
  // only required for the SphericalVertexCluster instantiation
  // TODO: Rename
  const double _parameter;

  // Input parameter
  const unsigned int _verticesPerCluster;
  const double       _relativeOverlap;
  const bool         _projectToInput;

  // Derived parameter
  double averageClusterRadius = 0;

  /// true if the mapping along some axis should be ignored
  /// has currently only dim x false entries, as integrated polynomials are irrelevant
  std::vector<bool> _deadAxis;

  /// polynomial treatment of the RBF system
  Polynomial _polynomial;

  /// @copydoc Mapping::mapConservative
  virtual void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  virtual void mapConsistent(DataID inputDataID, DataID outputDataID) override;

  /// export the center vertices of all clusters as a mesh with some additional data on it such as vertex count
  /// only enabled in debug builds and mainly for debugging purpose
  void exportClusterCentersAsVTU(mesh::Mesh &centers);
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::PartitionOfUnityMapping(
    Mapping::Constraint constraint,
    int                 dimension,
    double              parameter,
    Polynomial          polynomial,
    unsigned int        verticesPerCluster,
    double              relativeOverlap,
    bool                projectToInput)
    : Mapping(constraint, dimension),
      _parameter(parameter), _verticesPerCluster(verticesPerCluster), _relativeOverlap(relativeOverlap), _projectToInput(projectToInput), _polynomial(polynomial)
{
  PRECICE_ASSERT(this->getDimensions() <= 3);
  PRECICE_ASSERT(_relativeOverlap < 1, "The relative overlap has to be smaller than one.");
  PRECICE_ASSERT(_verticesPerCluster > 0, "The number of vertices per cluster has to be greater zero.");

  if (isScaledConsistent()) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
  // Simply set each axis to be active according to the space dimension.
  // We could even set each axis to be true, which would just lead to 'zero' columns in the polynomial QR.
  // Integrated polynomials are anyway not supported with PoU
  _deadAxis.clear();
  std::fill_n(std::back_inserter(_deadAxis), getDimensions(), false);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  PRECICE_TRACE();
  mesh::PtrMesh filterMesh, outMesh;
  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    filterMesh = this->output(); // remote
    outMesh    = this->input();  // local
  } else {
    filterMesh = this->input();  // remote
    outMesh    = this->output(); // local
  }

  if (outMesh->vertices().empty())
    return; // Ranks not at the interface should never hold interface vertices

  // TODO: Check again the tagging in combination with the partition construction (which mesh to use)
  // In order to construct the local partitions, we need all vertices with a distance of 2 x radius,
  // as the relevant partitions centers have a maximum distance of radius, and the proper construction of the
  // interpolant requires all vertices with a distance of radius from the center.
  auto bb = outMesh->getBoundingBox();

  if (averageClusterRadius == 0)
    averageClusterRadius = impl::estimatePartitionRadius(_verticesPerCluster, filterMesh, bb);

  // @TODO: This assert is not completely right, as it checks all dimensions for non-emptyness (which might not be the case).
  // However, with the current BB implementation, the expandBy function will just do nothing.
  PRECICE_ASSERT(!bb.empty());
  PRECICE_ASSERT(averageClusterRadius > 0);
  // Now we extend the bounding box by the radius
  bb.expandBy(1 * averageClusterRadius);

  // ... and tag all affected vertices
  auto verticesNew = filterMesh->index().getVerticesInsideBox(bb);

  std::for_each(verticesNew.begin(), verticesNew.end(), [&filterMesh](VertexID v) { filterMesh->vertices()[v].tag(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  // Probably nothing to be done here. There is no global ownership for matrix entries required and we tag all potentially locally relevant vertices already in the first round.
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  PRECICE_TRACE();

  precice::utils::Event e("map.pou.computeMapping.From" + this->input()->getName() + "To" + this->output()->getName(), precice::syncMode);

  // Recompute the whole clustering
  this->clear();

  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    inMesh  = this->output();
    outMesh = this->input();
  } else { // Consistent or scaled consistent
    inMesh  = this->input();
    outMesh = this->output();
  }

  // Step 1: get a tentative clustering consisting of centers and a radius from one of the available algorithms
  auto [averageClusterRadius_, centerCandidates] = impl::createUniformBlockPartitioning(inMesh, outMesh, _relativeOverlap, _verticesPerCluster, _projectToInput);

  averageClusterRadius = averageClusterRadius_;
  PRECICE_ASSERT(averageClusterRadius > 0 || inMesh->vertices().size() == 0 || outMesh->vertices().size() == 0);

  // Step 2: check, which of the resulting clusters are non-empty and register the cluster centers in a mesh
  // Here, the VertexCluster computes the matrix decompositions directly in case the cluster is non-empty
  mesh::Mesh centerMesh("clusterCentersMesh", this->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
  auto &     meshVertices = centerMesh.vertices();

  for (const auto &c : centerCandidates) {
    // We cannot simply take the vertex from the container, as the ID needs to match the cluster ID
    // That's required for the indexing and asserted below
    mesh::Vertex                                    center(c.getCoords(), meshVertices.size());
    SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T> p(inMesh->getDimensions(), center, averageClusterRadius, _parameter, _deadAxis, _polynomial, _verticesPerCluster, inMesh, outMesh);

    // Consider only non-empty clusters
    if (!p.isEmpty()) {
      PRECICE_ASSERT(center.getID() == _clusters.size(), center.getID(), _clusters.size());
      meshVertices.emplace_back(std::move(center));
      _clusters.emplace_back(std::move(p));
    }
  }
  // Log the average number of resulting clusters
  PRECICE_INFO("Number of total clusters (final): {}", _clusters.size());

  if (_clusters.size() > 0) {
    unsigned int nVertices   = std::accumulate(_clusters.begin(), _clusters.end(), static_cast<unsigned int>(0), [](auto &acc, auto &val) { return acc += val.getNumberOfInputVertices(); });
    unsigned int maxVertices = std::max_element(_clusters.begin(), _clusters.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices();
    unsigned int minVertices = std::min_element(_clusters.begin(), _clusters.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices();
    PRECICE_INFO("Average number of vertices per cluster {}", nVertices / _clusters.size());
    PRECICE_INFO("Maximum number of vertices per cluster {}", maxVertices);
    PRECICE_INFO("Minimum number of vertices per cluster {}", minVertices);
  }

  // Log a bounding box of the center mesh
  centerMesh.computeBoundingBox();
  PRECICE_INFO("Bounding Box of the cluster centers {}", centerMesh.getBoundingBox());

  // Step 3: index the clusters / the center mesh in order to create the output vertex -> cluster association
  query::Index clusterIndex(centerMesh);
  // Find all clusters the output vertex lies in, i.e., find all centers which have the distance of a cluster radius
  // Here: VertexID = clusterID
  // This could also be done on-the-fly in the map data phase, which would require to make the mesh as well as the indexTree member variables.
  PRECICE_DEBUG("Computing cluster-vertex association");
  for (const auto &vertex : outMesh->vertices()) {
    auto clusterIDs = clusterIndex.getVerticesInsideBox(vertex, averageClusterRadius);
    // Consider the case where we didn't find any partition (meshes don't match very well)
    if (clusterIDs.size() == 0) {
      PRECICE_WARN("Output vertex {} could not be assigned to a cluster. This means that the meshes probably do not match well geometry-wise.", vertex.getCoords());
      // TODO: Think about a proper way to handle this case, maybe set all radii to distance(v, closestvertex)?
      clusterIDs.emplace_back(clusterIndex.getClosestVertex(vertex.getCoords()).index);
    }

    // Step 4: compute the normalized weights of each output vertex for each partition
    PRECICE_ASSERT(clusterIDs.size() > 0, "No cluster found for vertex {}", vertex.getCoords());

    // Step 4a: compute the weight in each partition individually and store them in 'weights'
    std::vector<double> weights(clusterIDs.size());
    std::transform(clusterIDs.cbegin(), clusterIDs.cend(), weights.begin(), [&](const auto &ids) { return _clusters[ids].computeWeight(vertex); });
    double weightSum = std::accumulate(weights.begin(), weights.end(), static_cast<double>(0.));
    // TODO: This covers the edge case of vertices being at the edge of (several) clusters
    // In case the sum is equal to zero, we assign equal weights for all clusters
    if (!(weightSum > 0)) {
      PRECICE_ASSERT(weights.size() > 0);
      std::for_each(weights.begin(), weights.end(), [&weights](auto &w) { w = 1 / weights.size(); });
      weightSum = 1;
    }
    PRECICE_DEBUG("Weight sum {}", weightSum);
    PRECICE_DEBUG("Clusters {}", clusterIDs);
    PRECICE_DEBUG("V coords {}", vertex.getCoords());
    PRECICE_ASSERT(weightSum > 0);

    // Step 4b: scale the weight using the weight sum and store the normalized weight in all associated clusters
    for (unsigned int i = 0; i < clusterIDs.size(); ++i) {
      _clusters[clusterIDs[i]].setNormalizedWeight(weights[i] / weightSum, vertex.getID());
    }
  }

// Add a VTK export for visualization purposes
#ifndef NDEBUG
  exportClusterCentersAsVTU(centerMesh);
#endif

  // Set the computedMapping flag
  this->_hasComputedMapping = true;
  PRECICE_DEBUG("Compute Mapping is Completed.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _clusters.clear();
  // TODO: Don't reset this here
  averageClusterRadius      = 0;
  this->_hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConservative(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PRECICE_TRACE(inputDataID, outputDataID);
  // Execute the actual mapping evaluation in all clusters

  // 1. Reset the output data values as we need to accumulate data across clusters later on
  output()->data(outputDataID)->values().setZero();

  // 2. Iterate over all clusters and accumulate the result
  std::for_each(_clusters.begin(), _clusters.end(), [&](auto &p) { p.mapConservative(input()->data(inputDataID),
                                                                                     output()->data(outputDataID)); });
  // Set mapping finished
  std::for_each(_clusters.begin(), _clusters.end(), [&](auto &p) { p.setMappingFinished(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  // More detailed measurements
  PRECICE_TRACE(inputDataID, outputDataID);

  // 1. Reset the output data values as we need to accumulate data across clusters later on
  output()->data(outputDataID)->values().setZero();

  // 2. Execute the actual mapping evaluation in all vertex clusters and acccumulate the data
  std::for_each(_clusters.begin(), _clusters.end(), [&](auto &p) { p.mapConsistent(input()->data(inputDataID),
                                                                                   output()->data(outputDataID)); });

  // Set mapping finished
  std::for_each(_clusters.begin(), _clusters.end(), [&](auto &p) { p.setMappingFinished(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::exportClusterCentersAsVTU(mesh::Mesh &centerMesh)
{
  PRECICE_TRACE();

  auto dataRadius      = centerMesh.createData("radius", 1, -1);
  auto dataCardinality = centerMesh.createData("number-of-vertices", 1, -1);
  centerMesh.allocateDataValues();
  dataRadius->values().fill(averageClusterRadius);
  for (unsigned int i = 0; i < _clusters.size(); ++i) {
    dataCardinality->values()[i] = static_cast<double>(_clusters[i].getNumberOfInputVertices());
  }

  // We have to create the global offsets in order to export things in parallel
  if (utils::IntraComm::isSecondary()) {
    // send number of vertices
    PRECICE_DEBUG("Send number of vertices: {}", centerMesh.vertices().size());
    int numberOfVertices = centerMesh.vertices().size();
    utils::IntraComm::getCommunication()->send(numberOfVertices, 0);

    // receive vertex offsets
    mesh::Mesh::VertexOffsets vertexOffsets;
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets, 0);
    PRECICE_DEBUG("My vertex offsets: {}", vertexOffsets);
    PRECICE_ASSERT(centerMesh.getVertexOffsets().empty());
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  } else if (utils::IntraComm::isPrimary()) {

    mesh::Mesh::VertexOffsets vertexOffsets(utils::IntraComm::getSize());
    vertexOffsets[0] = centerMesh.vertices().size();

    // receive number of secondary vertices and fill vertex offsets
    for (int secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      int numberOfSecondaryRankVertices = -1;
      utils::IntraComm::getCommunication()->receive(numberOfSecondaryRankVertices, secondaryRank);
      PRECICE_ASSERT(numberOfSecondaryRankVertices >= 0);
      vertexOffsets[secondaryRank] = numberOfSecondaryRankVertices + vertexOffsets[secondaryRank - 1];
    }

    // broadcast vertex offsets
    PRECICE_DEBUG("My vertex offsets: {}", centerMesh.getVertexOffsets());
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets);
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  }

  io::ExportVTU exporter;
  exporter.doExport("pouCenters", "exports", centerMesh);
}
} // namespace mapping
} // namespace precice

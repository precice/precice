#pragma once

#include <Eigen/Core>
#include <numeric>

#include "com/Communication.hpp"
#include "io/ExportVTU.hpp"
#include "mapping/BatchedRBFSolver.hpp"
#include "mapping/impl/CreateClustering.hpp"
#include "mapping/impl/MappingDataCache.hpp"
#include "mapping/impl/SphericalVertexCluster.hpp"
#include "mesh/Filter.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "query/Index.hpp"
#include "utils/IntraComm.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

/**
 * Mapping using partition of unity decomposition strategies: The class here inherits from the Mapping
 * class and orchestrates the partitions (called vertex clusters) in order to represent a partition of unity.
 * This means in particular that the class computes the weights for the evaluation vertices and the necessary
 * association between evaluation vertices and the clusters during initialization and traverses through all
 * vertex clusters when evaluating the mapping.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class PartitionOfUnityMapping : public Mapping {
public:
  /**
   * Constructor, which mostly sets the mesh connectivity requirements and initializes member variables.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimension Dimensionality of the meshes
   * @param[in] function Radial basis function type used in interpolation
   * @param[in] polynomial The handling of the polynomial in the RBF system. Valid choices are 'off' and 'separate'
   * @param[in] verticesPerCluster Target number of vertices to be clustered together
   * @param[in] relativeOverlap Overlap between clusters: The parameter here determines the distance between two cluster
   * centers, given the cluster radius (already determined through \p verticesPerCluster ). A value of 1 would correspond
   * to no distance between cluster centers (i.e. completely overlapping clusters), 0 to distance of 2 x radius between
   * clusters centers.
   * @param[in] projectToInput if enabled, places the cluster centers at the closest vertex of the input mesh.
   * See also \ref mapping::impl::createClustering()
   */
  PartitionOfUnityMapping(
      Mapping::Constraint                   constraint,
      int                                   dimension,
      RADIAL_BASIS_FUNCTION_T               function,
      Polynomial                            polynomial,
      unsigned int                          verticesPerCluster,
      double                                relativeOverlap,
      bool                                  projectToInput,
      MappingConfiguration::GinkgoParameter ginkgoParameter          = MappingConfiguration::GinkgoParameter(),
      bool                                  computeEvaluationOffline = false);

  /**
   * Computes the clustering for the partition of unity method and fills the \p _clusters vector,
   * which allows to travers through all vertex cluster computed. Each vertex cluster in the vector
   * directly computes local mapping matrices and matrix decompositions.
   * In addition, the method computes the normalized weights (Shepard's method) for the partition
   * of unity method and stores them directly in each relevant vertex cluster.
   * In debug mode, the function also exports the partition centers as a separate mesh for visualization
   * purpose.
   */
  void computeMapping() final override;

  /// Clears a computed mapping by deleting the content of the \p _clusters vector.
  void clear() final override;

  /// tag the vertices required for the mapping
  void tagMeshFirstRound() final override;

  /// nothing to do here
  void tagMeshSecondRound() final override;

  /// name of the pum mapping
  std::string getName() const final override;

  void mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values) final override;

  /// the target values here remain unused, as we store the (intermediate) result directly in the cache
  void mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const Eigen::Ref<const Eigen::MatrixXd> &source, impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> target) final override;

  void updateMappingDataCache(impl::MappingDataCache &cache, const Eigen::Ref<const Eigen::VectorXd> &in) final override;

  void completeJustInTimeMapping(impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> buffer) final override;

  void initializeMappingDataCache(impl::MappingDataCache &cache) final override;

private:
  /// logger, as usual
  precice::logging::Logger _log{"mapping::PartitionOfUnityMapping"};

  void _computeCPU(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh, double clusterRadius, const std::vector<mesh::Vertex> &centerCandidates);

  /// main data container storing all the clusters, which need to be solved individually
  std::vector<SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>> _clusters;

  /// Radial basis function type used in interpolation
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// Input parameters provided by the user for the clustering algorithm:

  /// target number of input vertices for each cluster
  const unsigned int _verticesPerCluster;

  /// overlap of vertex clusters
  const double _relativeOverlap;

  /// toggles whether we project the cluster centers to the input mesh
  const bool _projectToInput;

  /// derived parameter based on the input above: the radius of each cluster
  double _clusterRadius = 0;

  /// polynomial treatment of the RBF system
  Polynomial _polynomial;

  const bool                                                 _useBatchedSolver;
  std::unique_ptr<BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>> _batchedSolver;

  MappingConfiguration::GinkgoParameter _ginkgoParameter;

  std::unique_ptr<mesh::Mesh> _centerMesh;

  /// @copydoc Mapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) override;

  /// Given an output vertex, computes the normalized PU weight at the given location
  /// The mesh name is only required for logging purposes
  /// returns the clusterIDs and all weights for these clusters
  std::pair<std::vector<int>, std::vector<double>> computeNormalizedWeight(const mesh::Vertex &v, std::string_view mesh);

  /// export the center vertices of all clusters as a mesh with some additional data on it such as vertex count
  /// only enabled in debug builds and mainly for debugging purpose
  void exportClusterCentersAsVTU(mesh::Mesh &centers);

  // Currently only valid for the Batched RBF solver, maybe move it to dedicated config struct
  const bool _computeEvaluationOffline;
};

template <typename RADIAL_BASIS_FUNCTION_T>
PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::PartitionOfUnityMapping(
    Mapping::Constraint                   constraint,
    int                                   dimension,
    RADIAL_BASIS_FUNCTION_T               function,
    Polynomial                            polynomial,
    unsigned int                          verticesPerCluster,
    double                                relativeOverlap,
    bool                                  projectToInput,
    MappingConfiguration::GinkgoParameter ginkgoParameter,
    bool                                  computeEvaluationOffline)
    : Mapping(constraint, dimension, false, Mapping::InitialGuessRequirement::None),
      _basisFunction(function), _verticesPerCluster(verticesPerCluster), _relativeOverlap(relativeOverlap),
      _projectToInput(projectToInput), _polynomial(polynomial), _useBatchedSolver(ginkgoParameter.executor != "cpu"),
      _ginkgoParameter(ginkgoParameter), _computeEvaluationOffline(computeEvaluationOffline)
{
  PRECICE_ASSERT(this->getDimensions() <= 3);
  PRECICE_ASSERT(_polynomial != Polynomial::ON, "Integrated polynomial is not supported for partition of unity data mappings.");
  PRECICE_ASSERT(_relativeOverlap < 1, "The relative overlap has to be smaller than one.");
  PRECICE_ASSERT(_verticesPerCluster > 0, "The number of vertices per cluster has to be greater zero.");
#ifdef PRECICE_NO_KOKKOS_KERNELS
  PRECICE_ASSERT(_useBatchedSolver == false, "Not implemented");
#endif
  if (isScaledConsistent()) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  PRECICE_TRACE();

  precice::profiling::Event e("map.pou.computeMapping.From" + this->input()->getName() + "To" + this->output()->getName(), profiling::Synchronize);

  /**
   * PARTITION OF UNITY (PUM) ALGORITHM OVERVIEW:
   *
   * The Partition of Unity Method decomposes a complex interpolation problem into multiple
   * overlapping local problems, each solved independently over smaller regions called clusters.
   * This makes the problem more manageable and improves computational efficiency.
   *
   * The algorithm consists of:
   * 1. Spatial clustering: Divide the domain into overlapping spherical clusters
   * 2. Local RBF systems: Each cluster builds its own RBF interpolant
   * 3. Weight computation: Assign normalized weights (Shepard's method) to blend local solutions
   * 4. Data mapping: Combine weighted contributions from all clusters containing a point
   */

  // Recompute the whole clustering
  PRECICE_ASSERT(!this->_hasComputedMapping, "Please clear the mapping before recomputing.");

  // Determine mesh roles based on mapping constraint type.
  // For conservative mappings: inMesh receives data, outMesh sends data (reversed roles).
  // For consistent mappings: inMesh is source, outMesh is target (normal roles).
  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    inMesh  = this->output();
    outMesh = this->input();
  } else { // Consistent or scaled consistent
    inMesh  = this->input();
    outMesh = this->output();
  }

  precice::profiling::Event eClusters("map.pou.computeMapping.createClustering.From" + this->input()->getName() + "To" + this->output()->getName());

  /**
   * STEP 1: SPATIAL CLUSTERING
   *
   * Create a set of overlapping spherical clusters that will partition the domain.
   * - Cluster centers are placed to optimally distribute vertices from the input mesh
   * - Cluster radius is computed to contain approximately _verticesPerCluster vertices
   * - Relative overlap controls cluster overlap: 0 (non-overlapping) to 1 (maximal overlap)
   *
   * Mathematical basis:
   * - Cluster radius r is determined such that each cluster contains ~verticesPerCluster points
   * - Distance between cluster centers is adjusted as: d = 2*r*(1-relativeOverlap)
   *   where relativeOverlap ∈ [0,1)
   */
  auto [clusterRadius, centerCandidates] = impl::createClustering(inMesh, outMesh, _relativeOverlap, _verticesPerCluster, _projectToInput);
  eClusters.stop();
  // Due to the aggressive filtering in the createClustering stage, the canddidates here are already final
  e.addData("n clusters", centerCandidates.size());

  _clusterRadius = clusterRadius;
  PRECICE_ASSERT(_clusterRadius > 0 || inMesh->nVertices() == 0 || outMesh->nVertices() == 0);

  if (_useBatchedSolver) {
    PRECICE_CHECK(!(outMesh->isJustInTime() || inMesh->isJustInTime()), "Just-in-time mappings are not implemented for Kokkos- or Ginkgo-based solvers.");
    precice::profiling::Event eBatched("map.pou.computeMapping.batchedSolver");
    _batchedSolver = std::make_unique<BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>>(_basisFunction, inMesh, outMesh,
                                                                                 centerCandidates, _clusterRadius,
                                                                                 _polynomial, _computeEvaluationOffline, _ginkgoParameter);

    // For the batched solver, we don't register the _centerMesh as such
    PRECICE_ASSERT(!_centerMesh, "The centerMesh is only utilized for the CPU variant");
  } else {
    _computeCPU(inMesh, outMesh, _clusterRadius, centerCandidates);
  }
  this->_hasComputedMapping = true;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::_computeCPU(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh, double clusterRadius, const std::vector<mesh::Vertex> &centerCandidates)
{
  /**
   * CPU VARIANT OF PUM COMPUTATION
   *
   * This method handles the single-threaded CPU execution of the partition of unity setup.
   * It creates SphericalVertexCluster objects that each:
   * - Build a local RBF interpolation matrix
   * - Compute matrix decompositions (Cholesky or QR) for efficient solving
   * - Store the decomposition for runtime evaluation
   */

  // STEP 1: CREATE CLUSTER INFRASTRUCTURE
  // Build a mesh to track cluster centers. Each cluster center becomes a mesh vertex.
  // Cluster ID = vertex ID in the center mesh (critical for indexing during mapping)
  _centerMesh        = std::make_unique<mesh::Mesh>("pou-centers-" + inMesh->getName(), this->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
  auto &meshVertices = _centerMesh->vertices();

  meshVertices.clear();
  _clusters.clear();
  _clusters.reserve(centerCandidates.size());

  // Iterate through all proposed cluster centers and create clusters.
  // CRITICAL INVARIANT: cluster ID must equal its index in _clusters vector.
  // This ensures O(1) lookup during the mapping phase later.
  for (const auto &c : centerCandidates) {
    // Create a new vertex for the cluster center with ID = current cluster count.
    // We cannot use c directly because its ID might not match our cluster vector index.
    const VertexID                                  vertexID = meshVertices.size();
    mesh::Vertex                                    center(c.getCoords(), vertexID);

    // Construct the SphericalVertexCluster which:
    // - Collects all input vertices within distance _clusterRadius
    // - Collects all output vertices within distance _clusterRadius
    // - Builds and decomposes the local RBF interpolation matrix
    SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T> cluster(center, _clusterRadius, _basisFunction, _polynomial, inMesh, outMesh);

    // Only keep non-empty clusters (a cluster with no vertices is unusable).
    // This safeguard handles degenerate geometry or extreme parameter combinations.
    if (!cluster.empty()) {
      PRECICE_ASSERT(center.getID() == static_cast<int>(_clusters.size()), center.getID(), _clusters.size());
      meshVertices.emplace_back(std::move(center));
      _clusters.emplace_back(std::move(cluster));
    }
  }

  // Log the average number of resulting clusters
  PRECICE_DEBUG("Partition of unity data mapping between mesh \"{}\" and mesh \"{}\": mesh \"{}\" on rank {} was decomposed into {} clusters.", this->input()->getName(), this->output()->getName(), inMesh->getName(), utils::IntraComm::getRank(), _clusters.size());

  if (_clusters.size() > 0) {
    PRECICE_DEBUG("Average number of vertices per cluster {}", std::transform_reduce(
                                                                   _clusters.begin(), _clusters.end(), size_t{0}, std::plus<>(),
                                                                   [](const auto &c) { return c.getNumberOfInputVertices(); }) /
                                                                   _clusters.size());
    PRECICE_DEBUG("Maximum number of vertices per cluster {}", std::max_element(_clusters.begin(), _clusters.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices());
    PRECICE_DEBUG("Minimum number of vertices per cluster {}", std::min_element(_clusters.begin(), _clusters.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices());
  }

  precice::profiling::Event eWeights("map.pou.computeMapping.computeWeights");
  // Log a bounding box of the center mesh
  _centerMesh->computeBoundingBox();
  PRECICE_DEBUG("Bounding Box of the cluster centers {}", _centerMesh->getBoundingBox());

  /**
   * STEP 2: COMPUTE PARTITION OF UNITY WEIGHTS (Shepard's Method)
   *
   * For each output vertex, we compute normalized weights for all clusters that contain it.
   * This is essential for blending the local RBF solutions from overlapping clusters.
   *
   * Mathematical formula (Shepard's method):
   * - For each cluster j, compute weight function w_j(x) = φ((R - dist(x, center_j))/R)
   *   where R = cluster radius, φ is a smooth cutoff (typically distance-based)
   * - Normalize: W_j(x) = w_j(x) / Σ_k w_k(x)
   * - This ensures Σ_j W_j(x) = 1 (partition of unity property)
   */
  PRECICE_DEBUG("Computing cluster-vertex association");
  for (const auto &vertex : outMesh->vertices()) {
    // For each output vertex, find all clusters containing it and compute weights.
    // This establishes the vertex-to-cluster mapping needed for runtime evaluation.
    auto [clusterIDs, normalizedWeights] = computeNormalizedWeight(vertex, outMesh->getName());

    // Store the normalized weight in each associated cluster.
    // Each cluster will use its weight to scale its local RBF contribution during mapping.
    for (unsigned int i = 0; i < clusterIDs.size(); ++i) {
      PRECICE_ASSERT(clusterIDs[i] < static_cast<int>(_clusters.size()));
      _clusters[clusterIDs[i]].setNormalizedWeight(normalizedWeights[i], vertex.getID());
    }
  }
  eWeights.stop();

  // Uncomment to add a VTK export of the cluster center distribution for visualization purposes
  // exportClusterCentersAsVTU(*_centerMesh);

  // we need the center mesh index data structure
  if (!outMesh->isJustInTime()) {
    _centerMesh.reset();
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::pair<std::vector<int>, std::vector<double>> PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::computeNormalizedWeight(const mesh::Vertex &vertex, std::string_view mesh)
{
  /**
   * COMPUTE NORMALIZED PARTITION OF UNITY WEIGHTS
   *
   * This function implements Shepard's interpolation weights for blending local RBF solutions.
   * For a query point (vertex), it:
   * 1. Finds all clusters whose domains contain the point
   * 2. Computes an individual weight for each cluster
   * 3. Normalizes weights to sum to 1 (partition of unity constraint)
   *
   * This enables smooth blending of overlapping local approximations.
   */

  // STEP 1: SPATIAL INDEXING
  // Use an R-tree spatial index built on cluster centers for efficient range queries.
  // The index structure allows O(log n) lookup instead of O(n) linear search.
  PRECICE_ASSERT(_centerMesh);
  query::Index &clusterIndex = _centerMesh->index();

  /**
   * STEP 2: RANGE QUERY FOR RELEVANT CLUSTERS
   *
   * Find all cluster centers that are exactly _clusterRadius away from the query vertex.
   * This defines which clusters' RBF interpolants contribute to the mapping at this point.
   *
   * By construction (see _computeCPU), cluster ID == index in _clusters vector,
   * so returned indices directly access the correct cluster objects.
   */
  auto       clusterIDs            = clusterIndex.getVerticesInsideBox(vertex, _clusterRadius);
  const auto localNumberOfClusters = clusterIDs.size();

  // Consider the case where we didn't find any cluster (meshes don't match very well)
  //
  // In principle, we could assign the vertex to the closest cluster using clusterIDs.emplace_back(clusterIndex.getClosestVertex(vertex.getCoords()).index);
  // However, this leads to a conflict with weights already set in the corresponding cluster, since we insert the ID and, later on, map the ID to a local weight index
  // Of course, we could rearrange the weights, but we want to avoid the case here anyway, i.e., prefer to abort.
  PRECICE_CHECK(localNumberOfClusters > 0,
                "Output vertex {} of mesh \"{}\" could not be assigned to any cluster in the rbf-pum mapping. This probably means that the meshes do not match well geometry-wise: Visualize the exported preCICE meshes to confirm. "
                "If the meshes are fine geometry-wise, you can try to increase the number of \"vertices-per-cluster\" (default is 50), the \"relative-overlap\" (default is 0.15), or disable the option \"project-to-input\". "
                "These options are only valid for the <mapping:rbf-pum-direct/> tag.",
                vertex.getCoords(), mesh);

  // Ensure we found at least one cluster (error condition if mesh doesn't match geometry)
  PRECICE_ASSERT(localNumberOfClusters > 0, "No cluster found for vertex {}", vertex.getCoords());

  /**
   * STEP 2b: COMPUTE LOCAL WEIGHTS
   *
   * For each relevant cluster, compute its individual weight using the cluster's weight function.
   * The weight function typically measures proximity to the cluster center with smooth falloff.
   * Common choice: w_j(x) = (R - dist(x, c_j))^3 (cubic function on [0, R])
   */
  std::vector<double> weights(localNumberOfClusters);
  std::transform(clusterIDs.cbegin(), clusterIDs.cend(), weights.begin(),
    [&](const auto &ids) { return _clusters[ids].computeWeight(vertex); });

  double weightSum = std::accumulate(weights.begin(), weights.end(), static_cast<double>(0.));

  /**
   * EDGE CASE HANDLING: Zero Weight Sum
   *
   * This occurs when a vertex lies exactly at cluster boundaries where all weights
   * simultaneously approach zero. For numerical stability, assign equal weights to all
   * contributing clusters. This maintains the partition of unity property.
   */
  if (weightSum <= 0) {
    PRECICE_ASSERT(weights.size() > 0);
    std::for_each(weights.begin(), weights.end(), [&weights](auto &w) { w = 1. / weights.size(); });
    weightSum = 1;
  }
  PRECICE_ASSERT(weightSum > 0);

  /**
   * STEP 2c: NORMALIZE WEIGHTS TO PARTITION OF UNITY
   *
   * Divide each weight by the sum: N_j(x) = w_j(x) / Σ_k w_k(x)
   * Result: Σ_j N_j(x) = 1 (partition of unity constraint guaranteed)
   *
   * This normalization ensures exact reproduction of constant functions and
   * enables proper blending of data at cluster boundaries.
   */
  std::transform(weights.begin(), weights.end(), weights.begin(),
    [weightSum](double w) { return w / weightSum; });

  return {clusterIDs, weights};
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  /**
   * CONSERVATIVE DATA MAPPING WITH PARTITION OF UNITY
   *
   * Maps data from source mesh to target mesh while conserving total quantity.
   * For each output vertex, aggregates weighted contributions from all overlapping clusters.
   *
   * Algorithm:
   * For each cluster c:
   *   - Evaluate RBF interpolant using cluster's input vertices
   *   - Weight result by normalized partition of unity weight
   *   - Accumulate weighted result into output vector
   *
   * Result: output_value = Σ_c W_c(x) * RBF_c(input_data)
   */
  PRECICE_TRACE();

  precice::profiling::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);

  // PRECONDITION: output must be zero before accumulation.
  // Each cluster contributes additively, so we need a clean slate.
  PRECICE_ASSERT(outData.isZero());

  if (_useBatchedSolver) {
    PRECICE_ASSERT(false, "Not implemented");
  } else {
    // Evaluate mapping by accumulating contributions from all clusters.
    // Each cluster:
    // 1. Evaluates its local RBF interpolant
    // 2. Applies its partition of unity weight
    // 3. Adds weighted result to output (accumulation)
    std::for_each(_clusters.begin(), _clusters.end(),
      [&](auto &cluster) { cluster.mapConservative(inData, outData); });
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  /**
   * CONSISTENT DATA MAPPING WITH PARTITION OF UNITY
   *
   * Maps data from source mesh to target mesh with pointwise interpolation.
   * Each output vertex receives data from the weighted blend of all cluster interpolants.
   *
   * Algorithm:
   * For each cluster c:
   *   - Evaluate RBF interpolant at output vertex locations
   *   - Weight result by normalized partition of unity weight
   *   - Accumulate weighted result into output vector
   *
   * Result: output_value(vertex) = Σ_c W_c(vertex) * RBF_c(input_data)
   *
   * Key difference from conservative:
   * - Consistent evaluates at output vertex positions (location-based)
   * - Conservative preserves global quantities (integral-based)
   */
  PRECICE_TRACE();

  precice::profiling::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);

  // PRECONDITION: output must be zero before accumulation.
  // Each cluster contributes additively to the final result.
  PRECICE_ASSERT(outData.isZero());

  if (_useBatchedSolver) {
    PRECICE_ASSERT(_batchedSolver, "Not initialized");
    _batchedSolver->solveConsistent(inData, outData);
  } else {
    // Evaluate mapping by accumulating contributions from all clusters.
    // Each cluster:
    // 1. Evaluates its local RBF interpolant at output vertex coordinates
    // 2. Applies its partition of unity weight
    // 3. Adds weighted result to output (accumulation)
    std::for_each(_clusters.begin(), _clusters.end(),
      [&](auto &clusters) { clusters.mapConsistent(inData, outData); });
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const Eigen::Ref<const Eigen::MatrixXd> &source,
                                                                         impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd>)
{
  precice::profiling::Event e("map.pou.mapConservativeAt.From" + input()->getName());
  // @todo: it would most probably be more efficient to first group the vertices we receive here according to the clusters and then compute the solution

  PRECICE_TRACE();
  PRECICE_ASSERT(_centerMesh);
  PRECICE_ASSERT(cache.p.size() == _clusters.size());
  PRECICE_ASSERT(cache.polynomialContributions.size() == _clusters.size());
  PRECICE_ASSERT(!_useBatchedSolver, "Not implemented"); // The PRECICE_CHECK is earlier

  mesh::Vertex vertex(coordinates.col(0), -1);
  for (Eigen::Index v = 0; v < coordinates.cols(); ++v) {
    vertex.setCoords(coordinates.col(v));
    auto [clusterIDs, normalizedWeights] = computeNormalizedWeight(vertex, this->input()->getName());
    // Use the weight to interpolate the solution
    for (std::size_t i = 0; i < clusterIDs.size(); ++i) {
      PRECICE_ASSERT(clusterIDs[i] < static_cast<int>(_clusters.size()));
      auto id = clusterIDs[i];
      // the input mesh refers here to a consistent constraint
      Eigen::VectorXd res = normalizedWeights[i] * source.col(v);
      _clusters[id].addWriteDataToCache(vertex, res, cache.polynomialContributions[id], cache.p[id], *this->output().get());
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::completeJustInTimeMapping(impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> buffer)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!cache.p.empty());
  PRECICE_ASSERT(!cache.polynomialContributions.empty());
  precice::profiling::Event e("map.pou.completeJustInTimeMapping.From" + input()->getName());

  for (std::size_t c = 0; c < _clusters.size(); ++c) {
    // If there is no contribution, we don't have to evaluate
    if (cache.p[c].squaredNorm() > 0) {
      _clusters[c].evaluateConservativeCache(cache.polynomialContributions[c], cache.p[c], buffer);
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::initializeMappingDataCache(impl::MappingDataCache &cache)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_hasComputedMapping);
  cache.p.resize(_clusters.size());
  cache.polynomialContributions.resize(_clusters.size());
  for (std::size_t c = 0; c < _clusters.size(); ++c) {
    _clusters[c].initializeCacheData(cache.polynomialContributions[c], cache.p[c], cache.getDataDimensions());
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::updateMappingDataCache(impl::MappingDataCache &cache, const Eigen::Ref<const Eigen::VectorXd> &in)
{
  // We cannot synchronize this event, as the call to this function is rank-local only
  precice::profiling::Event e("map.pou.updateMappingDataCache.From" + input()->getName());
  PRECICE_ASSERT(cache.p.size() == _clusters.size());
  PRECICE_ASSERT(cache.polynomialContributions.size() == _clusters.size());
  Eigen::Map<const Eigen::MatrixXd> inMatrix(in.data(), cache.getDataDimensions(), in.size() / cache.getDataDimensions());
  for (std::size_t c = 0; c < _clusters.size(); ++c) {
    _clusters[c].computeCacheData(inMatrix, cache.polynomialContributions[c], cache.p[c]);
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values)
{
  precice::profiling::Event e("map.pou.mapConsistentAt.From" + input()->getName());
  // @todo: it would most probably be more efficient to first group the vertices we receive here according to the clusters and then compute the solution
  PRECICE_TRACE();
  PRECICE_ASSERT(_centerMesh);
  PRECICE_ASSERT(!_useBatchedSolver, "Not implemented"); // The PRECICE_CHECK is earlier

  // First, make sure that everything is reset before we start
  values.setZero();

  mesh::Vertex vertex(coordinates.col(0), -1);
  for (Eigen::Index v = 0; v < values.cols(); ++v) {
    vertex.setCoords(coordinates.col(v));
    auto [clusterIDs, normalizedWeights] = computeNormalizedWeight(vertex, this->output()->getName());
    // Use the weight to interpolate the solution
    for (std::size_t i = 0; i < clusterIDs.size(); ++i) {
      PRECICE_ASSERT(clusterIDs[i] < static_cast<int>(_clusters.size()));
      auto id = clusterIDs[i];
      // the input mesh refers here to a consistent constraint
      Eigen::VectorXd localRes = normalizedWeights[i] * _clusters[id].interpolateAt(vertex, cache.polynomialContributions[id], cache.p[id], *this->input().get());
      values.col(v) += localRes;
    }
  }
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

  if (outMesh->empty())
    return; // Ranks not at the interface should never hold interface vertices

  // The geometric filter of the repartitioning is always disabled for the PU-RBF.
  // The main rationale: if we use only a fraction of the mesh then we might end up
  // with too few vertices per rank and we cannot prevent too few vertices from being
  // tagged, if we have filtered too much vertices beforehand. When using the filtering,
  // the user could increase the safety-factor or disable the filtering, but that's
  // a bit hard to understand for users. When no geometric filter is applid,
  // vertices().size() is here the same as getGlobalNumberOfVertices. Hence, it is much
  // safer to make use of the unfiltered mesh for the parallel tagging.
  //
  // Drawback: the "estimateClusterRadius" below makes use of the mesh R* index tree, and
  // constructing the tree on the (unfiltered) global mesh is computationally expensive (O( N logN)).
  // We could pre-filter the global mesh to a local fraction (using our own geometric filtering,
  // maybe with an increased safety margin or even an iterative increase of the safety margin),
  // but then there is again the question on how to do this in a safe way, without risking
  // failures depending on the partitioning. So we stick here to the computationally more
  // demanding, but safer version.
  // See also https://github.com/precice/precice/pull/1912#issuecomment-2551143620

  // Get the local bounding boxes
  auto localBB = outMesh->getBoundingBox();
  // we cannot check for empty'ness here, as a single output mesh vertex
  // would lead to a 0D box with zero volume (considered empty). Thus, we
  // simply check here for default'ness, which is equivalent to outMesh->empty()
  // further above
  PRECICE_ASSERT(!localBB.isDefault());

  if (_clusterRadius == 0)
    _clusterRadius = impl::estimateClusterRadius(_verticesPerCluster, filterMesh, localBB);

  PRECICE_DEBUG("Cluster radius estimate: {}", _clusterRadius);
  PRECICE_ASSERT(_clusterRadius > 0);

  // Now we extend the bounding box by the radius
  localBB.expandBy(2 * _clusterRadius);

  // ... and tag all affected vertices
  auto verticesNew = filterMesh->index().getVerticesInsideBox(localBB);

  std::for_each(verticesNew.begin(), verticesNew.end(), [&filterMesh](VertexID v) { filterMesh->vertex(v).tag(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  // Nothing to be done here. There is no global ownership for matrix entries required and we tag all potentially locally relevant vertices already in the first round.
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::exportClusterCentersAsVTU(mesh::Mesh &centerMesh)
{
  PRECICE_TRACE();

  auto dataRadius      = centerMesh.createData("radius", 1, -1);
  auto dataCardinality = centerMesh.createData("number-of-vertices", 1, -1);
  centerMesh.allocateDataValues();
  dataRadius->values().fill(_clusterRadius);
  for (unsigned int i = 0; i < _clusters.size(); ++i) {
    dataCardinality->values()[i] = static_cast<double>(_clusters[i].getNumberOfInputVertices());
  }

  // We have to create the global offsets in order to export things in parallel
  if (utils::IntraComm::isSecondary()) {
    // send number of vertices
    PRECICE_DEBUG("Send number of vertices: {}", centerMesh.nVertices());
    int numberOfVertices = centerMesh.nVertices();
    utils::IntraComm::getCommunication()->send(numberOfVertices, 0);

    // receive vertex offsets
    mesh::Mesh::VertexOffsets vertexOffsets;
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets, 0);
    PRECICE_DEBUG("Vertex offsets: {}", vertexOffsets);
    PRECICE_ASSERT(centerMesh.getVertexOffsets().empty());
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  } else if (utils::IntraComm::isPrimary()) {

    mesh::Mesh::VertexOffsets vertexOffsets(utils::IntraComm::getSize());
    vertexOffsets[0] = centerMesh.nVertices();

    // receive number of secondary vertices and fill vertex offsets
    for (int secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      int numberOfSecondaryRankVertices = -1;
      utils::IntraComm::getCommunication()->receive(numberOfSecondaryRankVertices, secondaryRank);
      PRECICE_ASSERT(numberOfSecondaryRankVertices >= 0);
      vertexOffsets[secondaryRank] = numberOfSecondaryRankVertices + vertexOffsets[secondaryRank - 1];
    }

    // broadcast vertex offsets
    PRECICE_DEBUG("Vertex offsets: {}", centerMesh.getVertexOffsets());
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets);
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  }

  dataRadius->setSampleAtTime(0, time::Sample{1, dataRadius->values()});
  dataCardinality->setSampleAtTime(0, time::Sample{1, dataCardinality->values()});
  io::ExportVTU exporter{"PoU", "exports", centerMesh, io::Export::ExportKind::TimeWindows, 1, utils::IntraComm::getRank(), utils::IntraComm::getSize()};
  exporter.doExport(0, 0.0);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _clusters.clear();
  // TODO: Don't reset this here
  _clusterRadius            = 0;
  this->_hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::string PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::getName() const
{
  return "partition-of-unity RBF";
}
} // namespace mapping
} // namespace precice

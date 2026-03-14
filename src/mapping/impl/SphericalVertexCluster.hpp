#pragma once

/**
 * @file SphericalVertexCluster.hpp
 * @brief A single spherical partition in the Partition-of-Unity (PUM) RBF mapping.
 *
 * ## What is a SphericalVertexCluster?
 * The Partition-of-Unity Method (PUM) decomposes a large interpolation problem
 * into many small, overlapping local sub-problems — each called a "cluster".
 * This class represents **one cluster**: a sphere of radius `_radius` centred at
 * `_center`, which:
 *   1. Collects all input mesh vertices inside the sphere  (_inputIDs).
 *   2. Collects all output mesh vertices inside the sphere (_outputIDs).
 *   3. Builds a dense local RBF interpolation system on those vertices.
 *   4. Stores Shepard (partition-of-unity) weights for the output vertices.
 *
 * ## Role in PartitionOfUnityMapping
 * `PartitionOfUnityMapping` owns a `std::vector<SphericalVertexCluster>` and
 * orchestrates them:
 *   - **computeMapping()** constructs each cluster and assigns PU weights.
 *   - **mapConsistent() / mapConservative()** iterate over all clusters and
 *     accumulate their weighted local RBF results into the global output.
 *
 * ## Data-Flow Summary
 *
 * ### Consistent mapping (interpolation):
 *   global_input[inputIDs] → local_in (nInput × nDims)
 *       → RBF solve → local_result (nOutput × nDims)
 *           → global_output[outputIDs] += result * normalizedWeights
 *
 * ### Conservative mapping (load-spreading):
 *   global_input[outputIDs] * normalizedWeights → local_in (nOutput × nDims)
 *       → RBF conservative solve → local_result (nInput × nDims)
 *           → global_output[inputIDs] += local_result
 *
 * ## Weighting Function
 * The per-cluster Shepard weight is evaluated using `CompactPolynomialC2`
 * (Wendland C2 kernel), giving smooth, non-negative weights that taper to zero
 * at the cluster boundary. Normalizing them across all clusters that cover a
 * given output vertex ensures the Partition-of-Unity (PU) property:
 *   Σ_c w_c(x) = 1  for every x in the domain.
 */

#include <Eigen/Core>

#include <boost/container/flat_set.hpp>

#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Filter.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"

namespace precice::mapping {

/**
 * The SphericalVertexCluster represents a single partition in the partition of unity mapping.
 * Hence, the PartitionOfUnity mapping class owns a vector of SphericalVertexClusters in order
 * to compute the mapping.
 * In its core, the class consists of a geometric center vertex and a radius representing the
 * spherical shape of the partition. In order to compute an RBF interpolant, the class stores
 * VertexIDs of the input mesh and the output mesh lying within the sphere and a
 * RadialBasisFctSolver to assemble and solve mapping matrices. The solver class here is exactly
 * the same class used in the plain RBF mapping.
 * Since each cluster maps data within its domain, the class is required to have similar
 * functions as the mapping classes in preCICE, i.e., mapConsistent and mapConservative for
 * the mapping execution as well as clear for the reset. These functions are always called
 * from the corresponding PartitionOfUnity mapping class, i.e.,
 * PartitionOfUnityMapping::mapConsistent calls the mapConsistent function of the (all elements
 * in the cluster vector) cluster here.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class SphericalVertexCluster {
public:
  /**
   * The constructor uses the index RTree of the input mesh and output mesh in order to collect
   * the vertexIDs of the input mesh and the output mesh lying within the spherical domain of the cluster.
   * Note that the index trees of the meshes are constructed in case they are empty.
   * If there are no input vertices or output vertices in the given domain ( \p center and \p radius ),
   * the cluster is considered empty ( see also \ref empty() ) and the constructor returns immediately.
   * If the cluster is non-empty, an RBF solver is constructed. The RBF solver assembles the mapping
   * matrices and computes the matrix decomposition directly.
   *
   * @param[in] center Spatial center of the vertex cluster
   * @param[in] radius Spatial radius of the cluster associated to the \p center
   * @param[in] function Radial basis function type used in interpolation
   * @param[in] polynomial The polynomial treatment in the RBF system.
   * @param[in] inputMesh mesh where the interpolants are build on, i.e., the input mesh for consistent
   *                      mappings and the output mesh for conservative mappings
   * @param[in] outputMesh mesh where we evaluate the interpolants, i.e., the output mesh consistent
   *                      mappings and the input mesh for conservative mappings
   */
  SphericalVertexCluster(mesh::Vertex            center,
                         double                  radius,
                         RADIAL_BASIS_FUNCTION_T function,
                         Polynomial              polynomial,
                         mesh::PtrMesh           inputMesh,
                         mesh::PtrMesh           outputMesh);

  /// Evaluates a conservative mapping and agglomerates the result in the given output data
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) const;

  /// Computes and saves the RBF coefficients
  void computeCacheData(const Eigen::Ref<const Eigen::MatrixXd> &globalIn, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coefficientsOut) const;

  /// Evaluates a consistent mapping and agglomerates the result in the given output data
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) const;

  /// Set the normalized weight for the given \p vertexID in the outputMesh
  void setNormalizedWeight(double normalizedWeight, VertexID vertexID);

  void evaluateConservativeCache(Eigen::MatrixXd &epsilon, const Eigen::MatrixXd &Au, Eigen::Ref<Eigen::MatrixXd> out);

  /// Compute the weight for a given vertex
  double computeWeight(const mesh::Vertex &v) const;

  Eigen::VectorXd interpolateAt(const mesh::Vertex &v, const Eigen::MatrixXd &poly, const Eigen::MatrixXd &coeffs, const mesh::Mesh &inMesh) const;
  /// Number of input vertices this partition operates on
  unsigned int getNumberOfInputVertices() const;

  /// The center coordinate of this cluster
  std::array<double, 3> getCenterCoords() const;

  /// Invalidates and erases data structures the cluster holds
  void clear();

  /// Returns, whether the current cluster is empty or not, where empty means that there
  /// are either no input vertices or output vertices.
  bool empty() const;

  void addWriteDataToCache(const mesh::Vertex &v, const Eigen::VectorXd &load, Eigen::MatrixXd &epsilon, Eigen::MatrixXd &Au,
                           const mesh::Mesh &inMesh);

  void initializeCacheData(Eigen::MatrixXd &polynomial, Eigen::MatrixXd &coeffs, const int nComponents);

private:
  /// logger, as usual
  mutable precice::logging::Logger _log{"mapping::SphericalVertexCluster"};

  /// center vertex of the cluster
  mesh::Vertex _center;

  /// radius of the vertex cluster
  const double _radius;

  /// The RBF solver
  RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> _rbfSolver;

  // Stores the global IDs of the vertices so that we can apply a binary
  // search in order to query specific objects. Here we have logarithmic
  // complexity for setting the weights (only performed once), and constant
  // complexity when traversing through the IDs (performed in each iteration).

  /// Global VertexIDs associated to the input mesh (see constructor)
  boost::container::flat_set<VertexID> _inputIDs;

  /// Global VertexIDs associated to the output mesh (see constructor)
  boost::container::flat_set<VertexID> _outputIDs;

  /// Vector containing the normalized weights used for the output mesh data
  /// (consistent mapping) or input mesh data (conservative data)
  Eigen::VectorXd _normalizedWeights;

  /// Polynomial treatment in the RBF solver
  Polynomial _polynomial;

  RADIAL_BASIS_FUNCTION_T _function;

  /// The weighting function
  CompactPolynomialC2 _weightingFunction;

  /// Boolean switch in order to indicate that a mapping was computed
  bool _hasComputedMapping = false;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::SphericalVertexCluster(
    mesh::Vertex            center,
    double                  radius,
    RADIAL_BASIS_FUNCTION_T function,
    Polynomial              polynomial,
    mesh::PtrMesh           inputMesh,
    mesh::PtrMesh           outputMesh)
    : _center(center), _radius(radius), _polynomial(polynomial), _function(function), _weightingFunction(radius)
{
  PRECICE_TRACE(_center.getCoords(), _radius);
  precice::profiling::Event eq("map.pou.computeMapping.queryVertices");
  // Disable integrated polynomial, as it might cause locally singular matrices
  PRECICE_ASSERT(_polynomial != Polynomial::ON, "Integrated polynomial is not supported for partition of unity data mappings.");

  // Get vertices to be mapped
  // Subtract a safety margin to exclude the vertices at the edge
  auto outIDs = outputMesh->index().getVerticesInsideBox(center, radius - math::NUMERICAL_ZERO_DIFFERENCE);
  // Constructing the partition when we don't have evaluation points is pointless
  auto inIDs = inputMesh->index().getVerticesInsideBox(center, radius);

  // Transform the vector to the appropriate boost data structure
  // The IDs are sorted in the boost flat_set, hence, the function here has N log(N) complexity
  _inputIDs.insert(inIDs.begin(), inIDs.end());
  _outputIDs.insert(outIDs.begin(), outIDs.end());
  eq.stop();
  // If the cluster is empty, we return immediately
  if (empty()) {
    return;
  }

  PRECICE_DEBUG("SphericalVertexCluster input size: {}", inIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster output size: {}", outIDs.size());

  // The polynomial system is underdetermined if inIDs.size() < dimension + 1. However, the dynamic adoption of the axis in the RBF solver
  // disables axis, if necessary. Hence, we don't disable the complete polynomial here for underdetermined systems. The case should anyway
  // only occur for an almost unreasonable small vertices-per-cluster configuration.

  // Construct the solver. Here, the constructor of the RadialBasisFctSolver computes already the decompositions etc, such that we can mark the
  // mapping in this cluster as computed (mostly for debugging purpose)
  std::vector<bool>         deadAxis(inputMesh->getDimensions(), false);
  precice::profiling::Event e("map.pou.computeMapping.rbfSolver");
  _rbfSolver          = RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{function, *inputMesh.get(), _inputIDs, *outputMesh.get(), _outputIDs, deadAxis, _polynomial};
  _hasComputedMapping = true;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::setNormalizedWeight(double normalizedWeight, VertexID id)
{
  PRECICE_ASSERT(_outputIDs.size() > 0);
  // The given global vertex ID must belong to this cluster's output set.
  PRECICE_ASSERT(_outputIDs.contains(id), id);
  // Normalized weights must be strictly positive — zero would mean the cluster
  // contributes nothing at this vertex (it should not be in this cluster at all).
  PRECICE_ASSERT(normalizedWeight > 0);

  // Lazy initialization: allocate the weight vector on first use.
  // Size = number of output vertices in this cluster (local indexing, not global).
  if (_normalizedWeights.size() == 0)
    _normalizedWeights.resize(_outputIDs.size());

  // boost::container::flat_set stores elements in sorted order and provides
  // O(log N) lookup. `index_of` converts the iterator to a 0-based position,
  // giving us the local index that matches `_normalizedWeights` ordering.
  // This is the only time we pay O(log N) per vertex — the actual mapping
  // loop uses sequential access (O(1) per step via nth()).
  auto localID = _outputIDs.index_of(_outputIDs.find(id));

  PRECICE_ASSERT(static_cast<Eigen::Index>(localID) < _normalizedWeights.size(), localID, _normalizedWeights.size());
  _normalizedWeights[localID] = normalizedWeight;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) const
{
  /**
   * LOCAL CONSERVATIVE MAPPING FOR ONE CLUSTER
   *
   * Conservative mapping means the *total quantity* (integral / sum) is
   * preserved across the mapping. Mathematically:
   *   Σ_j outData[j] == Σ_i inData[i]  (weighted by the PU weights)
   *
   * The cluster's contribution to the global conservative solve is:
   *   out[inputIDs] += C^{-T} * A^T * (inData[outputIDs] * w)   (RBF formula)
   * where w = normalizedWeights (Shepard PU weights for this cluster).
   *
   * Key: weights are applied to the INPUT data (output mesh side) BEFORE
   * the solve, so each cluster only "sees" its share of the total load.
   * Multiple clusters accumulate into the global outData with +=.
   */

  // Guards: empty partitions are never stored; mapping must be pre-computed.
  PRECICE_ASSERT(!empty());
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(_normalizedWeights.size() == static_cast<Eigen::Index>(_outputIDs.size()));

  const unsigned int nComponents = inData.dataDims; // e.g. 1 for scalar, 3 for 3D vector
  const auto        &localInData = inData.values;   // flat global array: vertex * nComponents

  // Allocate a dense local matrix (nOutputInCluster × nComponents).
  // `getOutputSize()` == _outputIDs.size() (+ polynomial rows, already zero-padded).
  // TODO: We can probably reduce the temporary allocations here
  Eigen::MatrixXd in(_rbfSolver.getOutputSize(), nComponents);

  // STEP 1: Extract & weight the relevant portion of the global input data.
  // For conservative mapping, the "input" to the local solver is the OUTPUT mesh data
  // (note the reversed naming: the global outMesh == local inputMesh for conservative).
  // We multiply each row by the normalized Shepard weight so that overlapping clusters
  // collectively sum to the original unweighted value.
  for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
    for (unsigned int c = 0; c < nComponents; ++c) {
      const auto dataIndex = *(_outputIDs.nth(i)); // global vertex ID → flat array offset
      PRECICE_ASSERT(dataIndex * nComponents + c < localInData.size(), dataIndex * nComponents + c, localInData.size());
      PRECICE_ASSERT(_normalizedWeights[i] > 0, _normalizedWeights[i], i);
      // Apply partition-of-unity weight: each cluster only contributes its
      // fraction of the load. This ensures Σ_clusters contribution = full load.
      in(i, c) = localInData[dataIndex * nComponents + c] * _normalizedWeights[i];
    }
  }

  // STEP 2: Solve the local RBF system conservatively.
  // Internally: result = C^{-1} * A^T * in  (conservative solve via RadialBasisFctSolver)
  // result has shape (nInputInCluster × nComponents).
  Eigen::MatrixXd result = _rbfSolver.solveConservative(in, _polynomial);
  PRECICE_ASSERT(result.rows() == static_cast<Eigen::Index>(_inputIDs.size()));

  // STEP 3: Accumulate the local result into the global output vector.
  // += is critical: multiple clusters will add their contributions to the same
  // global input vertex, and their sum gives the final conservative mapped value.
  for (unsigned int i = 0; i < _inputIDs.size(); ++i) {
    for (unsigned int c = 0; c < nComponents; ++c) {
      const auto dataIndex = *(_inputIDs.nth(i)); // local row i → global vertex ID
      PRECICE_ASSERT(dataIndex * nComponents + c < outData.size(), dataIndex * nComponents + c, outData.size());
      outData[dataIndex * nComponents + c] += result(i, c); // accumulate
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::evaluateConservativeCache(Eigen::MatrixXd &epsilon, const Eigen::MatrixXd &Au, Eigen::Ref<Eigen::MatrixXd> out)
{
  Eigen::MatrixXd localIn(_inputIDs.size(), Au.cols());
  _rbfSolver.evaluateConservativeCache(epsilon, Au, localIn);
  // Step 3: now accumulate the result into our global output data
  for (std::size_t i = 0; i < _inputIDs.size(); ++i) {
    const auto dataIndex = *(_inputIDs.nth(i));
    PRECICE_ASSERT(dataIndex < out.cols(), out.cols());
    out.col(dataIndex) += localIn.row(i);
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::computeCacheData(const Eigen::Ref<const Eigen::MatrixXd> &globalIn, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coeffOut) const
{
  PRECICE_TRACE();
  Eigen::MatrixXd in(_rbfSolver.getInputSize(), globalIn.rows());
  // Step 1: extract the relevant input data from the global input data and store
  // it in a contiguous array, which is required for the RBF solver (last polyparams entries remain zero)
  for (std::size_t i = 0; i < _inputIDs.size(); i++) {
    const auto dataIndex = *(_inputIDs.nth(i));
    PRECICE_ASSERT(dataIndex < globalIn.cols(), globalIn.cols());
    in.row(i) = globalIn.col(dataIndex);
  }
  // Step 2: solve the system using a consistent constraint
  _rbfSolver.computeCacheData(in, _polynomial, polyOut, coeffOut);
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::interpolateAt(const mesh::Vertex &v, const Eigen::MatrixXd &poly, const Eigen::MatrixXd &coeffs, const mesh::Mesh &inMesh) const
{
  PRECICE_TRACE();
  return _rbfSolver.interpolateAt(v, poly, coeffs, _function, _inputIDs, inMesh);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::addWriteDataToCache(const mesh::Vertex &v, const Eigen::VectorXd &load,
                                                                          Eigen::MatrixXd &epsilon, Eigen::MatrixXd &Au,
                                                                          const mesh::Mesh &inMesh)
{
  PRECICE_TRACE();
  _rbfSolver.addWriteDataToCache(v, load, epsilon, Au, _function, _inputIDs, inMesh);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::initializeCacheData(Eigen::MatrixXd &poly, Eigen::MatrixXd &coeffs, const int nComponents)
{
  /**
   * INITIALIZE JIT CACHE BUFFERS FOR THIS CLUSTER
   *
   * Called once before the just-in-time (JIT) mapping loop begins.
   * Pre-allocates the two buffers that will hold the precomputed RBF solution
   * for one time step, so that per-query evaluations (interpolateAt /
   * addWriteDataToCache) can proceed without further allocation.
   *
   * - `poly`  : polynomial coefficients matrix.
   *             Only meaningful when polynomial=SEPARATE; otherwise left at
   *             size 0 (checked in the solver via PRECICE_ASSERT).
   *             Shape: (nPolynomialParams × nComponents)
   *
   * - `coeffs`: RBF lambda coefficients ("C^{-1} * inputData" for consistent).
   *             Shape: (nInputVerticesInCluster × nComponents)
   *             Each column corresponds to one data component.
   */
  if (Polynomial::SEPARATE == _polynomial) {
    // Allocate polynomial coefficient buffer: one row per polynomial basis term
    // (1 constant + number_of_active_dimensions linear terms).
    poly.resize(_rbfSolver.getNumberOfPolynomials(), nComponents);
  }
  // Allocate RBF coefficient buffer: one row per input vertex in this cluster.
  coeffs.resize(_inputIDs.size(), nComponents);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) const
{
  /**
   * LOCAL CONSISTENT MAPPING FOR ONE CLUSTER
   *
   * Consistent (interpolation) mapping evaluates the RBF interpolant at output
   * vertex positions. For this cluster, the formula is:
   *
   *   outData[outputIDs] += w * (A * C^{-1} * inData[inputIDs])
   *
   * where:
   *   - C = RBF kernel matrix at input vertices (assembled in constructor)
   *   - A = evaluation matrix at output vertices (assembled in constructor)
   *   - w = normalizedWeights (Shepard PU weight for this cluster)
   *
   * The weight `w` is applied AFTER interpolation (to the result), so each
   * cluster contributes a fraction of the final interpolated value. Summing
   * contributions from all clusters that cover an output vertex gives the
   * final mapped value (guaranteed by the PU normalization: Σ w_c = 1).
   *
   * Multiple clusters accumulate with += into the global outData.
   */

  // Guards: empty partitions are never stored; mapping must be pre-computed.
  PRECICE_ASSERT(!empty());
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(_normalizedWeights.size() == static_cast<Eigen::Index>(_outputIDs.size()));

  const unsigned int nComponents = inData.dataDims; // e.g. 1 for scalar, 3 for 3D vector
  const auto        &localInData = inData.values;   // flat global array: vertex * nComponents

  // Allocate the local input matrix (nInputInCluster × nComponents).
  // getInputSize() >= _inputIDs.size() because the solver may have extra rows
  // for the polynomial system (those rows stay zero for the RBF-only part).
  Eigen::MatrixXd in(_rbfSolver.getInputSize(), nComponents);

  // STEP 1: Extract the relevant global input data into a dense local matrix.
  // We only pull values for vertices that belong to this cluster (_inputIDs),
  // ignoring all other global vertices. The polynomial rows (beyond inputIDs)
  // remain zero-initialized — the solver expects zeros there.
  for (unsigned int i = 0; i < _inputIDs.size(); i++) {
    for (unsigned int c = 0; c < nComponents; ++c) {
      const auto dataIndex = *(_inputIDs.nth(i)); // global vertex ID → flat array offset
      PRECICE_ASSERT(dataIndex * nComponents + c < localInData.size(), dataIndex * nComponents + c, localInData.size());
      in(i, c) = localInData[dataIndex * nComponents + c];
    }
  }

  // STEP 2: Solve the local RBF system (consistent).
  // Internally: coeffs = C^{-1} * in, result = A * coeffs
  // result has shape (nOutputInCluster × nComponents).
  Eigen::MatrixXd result = _rbfSolver.solveConsistent(in, _polynomial);
  PRECICE_ASSERT(static_cast<Eigen::Index>(_outputIDs.size() * nComponents) == result.size());

  // STEP 3: Scale each local result by this cluster's PU weight, then accumulate
  // into the global output vector at the correct global vertex positions.
  // The PU weight ensures that overlapping cluster contributions sum to 1.
  for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
    for (unsigned int c = 0; c < nComponents; ++c) {
      const auto dataIndex = *(_outputIDs.nth(i)); // local row i → global vertex ID
      PRECICE_ASSERT(dataIndex * nComponents + c < outData.size(), dataIndex * nComponents + c, outData.size());
      PRECICE_ASSERT(_normalizedWeights[i] > 0);
      // Multiply by the normalized Shepard weight before accumulation.
      // When all clusters have done this, the sum at each output vertex equals
      // the fully weighted interpolated value (PU property guaranteed).
      outData[dataIndex * nComponents + c] += result(i, c) * _normalizedWeights[i];
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
double SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::computeWeight(const mesh::Vertex &v) const
{
  /**
   * SHEPARD WEIGHT COMPUTATION
   *
   * Returns the unnormalized weight w_c(v) for vertex v in this cluster.
   * The weight is evaluated using the cluster's `_weightingFunction`, which
   * is a `CompactPolynomialC2` (Wendland C2) kernel with support radius = _radius:
   *
   *   w_c(v) = phi_C2(dist(v, center) / radius)
   *          = (1 - r)^4 * (4r + 1),   r = dist / radius
   *
   * Properties:
   * - w_c(v) == 0 exactly when dist(v, center) >= radius  (compact support).
   * - w_c(v) > 0 for all v strictly inside the sphere.
   * - Smooth (C^2): two continuous derivatives, ensuring smooth blending.
   *
   * These weights are later normalized by PartitionOfUnityMapping::computeNormalizedWeight():
   *   W_c(v) = w_c(v) / Σ_{k} w_k(v),  so that Σ_c W_c(v) = 1 (PU property).
   *
   * Note: dead axes are not considered here — because all 3 components of the
   * 3D coordinates are used, and a dead axis has identical values for all vertices,
   * so the coordinate difference along that axis is 0 anyway.
   */
  // Compute squared Euclidean distance from query vertex to the cluster center.
  // We use all 3 components (dead axes contribute 0 difference automatically).
  auto res = computeSquaredDifference(_center.rawCoords(), v.rawCoords(), {{true, true, true}});
  // Evaluate the Wendland C2 weight function at this distance.
  return _weightingFunction.evaluate(std::sqrt(res));
}

template <typename RADIAL_BASIS_FUNCTION_T>
unsigned int SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::getNumberOfInputVertices() const
{
  return _inputIDs.size();
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::array<double, 3> SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::getCenterCoords() const
{
  return _center.rawCoords();
}

template <typename RADIAL_BASIS_FUNCTION_T>
bool SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::empty() const
{
  return _inputIDs.size() == 0;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _inputIDs.clear();
  _outputIDs.clear();
  _rbfSolver.clear();
  _hasComputedMapping = false;
}
} // namespace precice::mapping

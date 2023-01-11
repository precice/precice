#pragma once

#include <Eigen/Core>

#include <boost/container/flat_set.hpp>

#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"

namespace precice {

namespace mapping {

/**
 * The SphericalVertexCluster represents a single partition in the partition of unity mapping.
 * Hence, the PartitionOfUnity mapping class owns a vector of SphericalVertexClusters in order
 * to compute the mapping.
 * In its core, the class consists of a geometric center vertex and a radius representing the
 * spherical shape of the partition. In order to compute an RBF interpolant, the class stores
 * VertexIDs of the input mesh and the output mesh lying within the sphere and a
 * RadialBasisFctSolver to assemble and solve mapping matrices. The solver class here is exactly
 * the same class used in the plain RBF mapping.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class SphericalVertexCluster {
public:
  /**
   * The constructor uses the index RTree of the input mesh and output mesh in order to collect
   * the vertexIDs of the input mesh and the output mesh lying within the spherical domain of the cluster.
   * Note that the index trees of the meshes are constructed in case they are empty.
   * If there are no input vertices or output verices in the given domain ( \p center and \p radius ),
   * the cluster is considered empty ( see also \ref empty() ) and the constructor returns immediately.
   * If the cluster is non-empty, an RBF solver is constructed. The RBF solver assembles the mapping
   * matrices and computes the matrix decomposition directly.
   *
   * @param[in] center Spatial center of the vertex cluster
   * @param[in] radius Spatial radius of the cluster associated to the \p center
   * @param[in] parameter Shape parameter or support radius of the RBF (TODO)
   * @param[in] deadAxis dead axis as set by the user. Required for the RBF solver
   * @param[in] polynomial The polynomial treatment in the RBF system.
   * @param[in] inputMesh mesh where the interpolants are build on, i.e., the input mesh for consistent
   *                      mappings and the output mesh for conservative mappings
   * @param[in] outputMesh mesh where we evaluate the interpolants, i.e., the output mesh consistent
   *                      mappings and the input mesh for conservative mappings
   */
  SphericalVertexCluster(mesh::Vertex      center,
                         double            radius,
                         double            parameter,
                         std::vector<bool> deadAxis,
                         Polynomial        polynomial,
                         mesh::PtrMesh     inputMesh,
                         mesh::PtrMesh     outputMesh);

  /// Evaluates a conservative mapping and agglomerates the result in the given output data
  void mapConservative(mesh::PtrData inputData, mesh::PtrData outputData) const;

  /// Evaluates a consistent mapping and agglomerates the result in the given output data
  void mapConsistent(mesh::PtrData inputData, mesh::PtrData outputData) const;

  /// Set the normalized weight for the given \p vertexID in the outputMesh
  void setNormalizedWeight(double normalizedWeight, VertexID vertexID);

  /// Compute the weight for a given vertex
  double computeWeight(const mesh::Vertex &v) const;

  /// Number of input vertices this partition operates on
  unsigned int getNumberOfInputVertices() const;

  /// The center coordinate of this cluster
  std::array<double, 3> getCenterCoords() const;

  /// Invalidates and erases data structures the cluster holds
  void clear();

  /// Returns, whether the current cluster is empty or not, where empty means that there
  /// are either no input vertices or output vertices.
  bool empty() const;

private:
  /// logger, as usual
  precice::logging::Logger _log{"mapping::SphericalVertexCluster"};

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

  /// The weighting function
  CompactPolynomialC2 _weightingFunction;

  /// Boolean switch in order to indicate that a mapping was computed
  bool _hasComputedMapping = false;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::SphericalVertexCluster(
    mesh::Vertex      center,
    double            radius,
    double            parameter,
    std::vector<bool> deadAxis,
    Polynomial        polynomial,
    mesh::PtrMesh     inputMesh,
    mesh::PtrMesh     outputMesh)
    : _center(center), _radius(radius), _polynomial(polynomial), _weightingFunction(radius)
{
  // Disable integrated polynomial, as it might cause locally singular matrices
  PRECICE_ASSERT(_polynomial != Polynomial::ON, "Integrated polynomial is not supported for partition of unity data mappings.")
  PRECICE_ASSERT(deadAxis.size() == inputMesh->getDimensions());
  PRECICE_DEBUG("Center coordinates: {}", _center.getCoords());
  PRECICE_DEBUG("SphericalVertexCluster radius: {}", _radius);

  // Get vertices to be mapped
  auto outIDs = outputMesh->index().getVerticesInsideBox(center, radius);
  // Constructing the partition when we don't have evaluation points is pointless
  auto inIDs = inputMesh->index().getVerticesInsideBox(center, radius);

  // Transform the vector to the appropriate boost data structure
  // The IDs are sorted in the boost flat_set, hence, the function here has N log(N) complexity
  _inputIDs.insert(inIDs.begin(), inIDs.end());
  _outputIDs.insert(outIDs.begin(), outIDs.end());

  // If the cluster is empty, we return immediately
  if (empty()) {
    return;
  }

  PRECICE_DEBUG("SphericalVertexCluster input size: {}", inIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster output size: {}", outIDs.size());

  // TODO
  // Otherwise the system is underdetermined
  // if (inIDs.size() < dimension + 1)
  //   _polynomial = Polynomial::OFF;

  // Construct the solver
  _rbfSolver          = RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{RADIAL_BASIS_FUNCTION_T{parameter}, *inputMesh.get(), _inputIDs, *outputMesh.get(), _outputIDs, deadAxis, _polynomial};
  _hasComputedMapping = true;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::setNormalizedWeight(double normalizedWeight, VertexID id)
{
  PRECICE_ASSERT(_outputIDs.size() > 0);
  PRECICE_ASSERT(_outputIDs.contains(id), id);
  PRECICE_ASSERT(normalizedWeight > 0);

  if (_normalizedWeights.size() == 0)
    _normalizedWeights.resize(_outputIDs.size());

  // The find method of boost flat_set comes with O(log(N)) complexity (the more expensive part here)
  auto localID = _outputIDs.index_of(_outputIDs.find(id));

  PRECICE_ASSERT(localID < _normalizedWeights.size(), localID, _normalizedWeights.size());
  _normalizedWeights[localID] = normalizedWeight;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConservative(mesh::PtrData inputData, mesh::PtrData outputData) const
{
  // First, a few sanity checks. Empty partitions shouldn't be stored at all
  PRECICE_ASSERT(!empty());
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(_normalizedWeights.size() == _outputIDs.size());

  // Define an alias for data dimension in order to avoid ambiguity
  const unsigned int nComponents = inputData->getDimensions();
  const auto &       localInData = inputData->values();

  // TODO: We can probably reduce the temporary allocations here
  Eigen::VectorXd in(_rbfSolver.getEvaluationMatrix().rows());

  // Now we perform the data mapping component-wise
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Step 1: extract the relevant input data from the global input data and store
    // it in a contiguous array, which is required for the RBF solver
    for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
      const auto dataIndex = *(_outputIDs.nth(i));
      PRECICE_ASSERT(dataIndex * nComponents + c < localInData.size(), dataIndex * nComponents + c, localInData.size());
      PRECICE_ASSERT(_normalizedWeights[i] > 0, _normalizedWeights[i], i);
      // here, we also directly apply the weighting, i.e., we split the input data
      in[i] = localInData[dataIndex * nComponents + c] * _normalizedWeights[i];
    }

    // Step 2: solve the system using a conservative constraint
    auto result = _rbfSolver.solveConservative(in, _polynomial);
    PRECICE_ASSERT(result.size() == _inputIDs.size());

    // Step 3: now accumulate the result into our global output data
    for (unsigned int i = 0; i < _inputIDs.size(); ++i) {
      const auto dataIndex = *(_inputIDs.nth(i));
      PRECICE_ASSERT(dataIndex * nComponents + c < outputData->values().size(), dataIndex * nComponents + c, outputData->values().size());
      outputData->values()[dataIndex * nComponents + c] += result(i);
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConsistent(mesh::PtrData inputData, mesh::PtrData outputData) const
{
  // First, a few sanity checks. Empty partitions shouldn't be stored at all
  PRECICE_ASSERT(!empty());
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(_normalizedWeights.size() == _outputIDs.size());

  // Define an alias for data dimension in order to avoid ambiguity
  const unsigned int nComponents = inputData->getDimensions();
  const auto &       localInData = inputData->values();

  Eigen::VectorXd in(_rbfSolver.getEvaluationMatrix().cols());

  // Now we perform the data mapping component-wise
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Step 1: extract the relevant input data from the global input data and store
    // it in a contiguous array, which is required for the RBF solver (last polyparams entries remain zero)
    for (unsigned int i = 0; i < _inputIDs.size(); i++) {
      const auto dataIndex = *(_inputIDs.nth(i));
      PRECICE_ASSERT(dataIndex * nComponents + c < localInData.size(), dataIndex * nComponents + c, localInData.size());
      in[i] = localInData[dataIndex * nComponents + c];
    }

    // Step 2: solve the system using a consistent constraint
    auto result = _rbfSolver.solveConsistent(in, _polynomial);
    PRECICE_ASSERT(_outputIDs.size() == result.size());

    // Step 3: now accumulate the result into our global output data
    for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
      const auto dataIndex = *(_outputIDs.nth(i));
      PRECICE_ASSERT(dataIndex * nComponents + c < outputData->values().size(), dataIndex * nComponents + c, outputData->values().size());
      PRECICE_ASSERT(_normalizedWeights[i] > 0);
      // here, we also directly apply the weighting, i.e., split the result data
      outputData->values()[dataIndex * nComponents + c] += result(i) * _normalizedWeights[i];
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
double SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::computeWeight(const mesh::Vertex &v) const
{
  // Assume that the local interpolant and the weighting function are the same
  // TODO: We don't need to reduce the dead coordinates here as the values should reduce anyway
  auto res = computeSquaredDifference(_center.rawCoords(), v.rawCoords(), {{true, true, true}});
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
  return _outputIDs.size() == 0 || _inputIDs.size() == 0;
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
} // namespace mapping
} // namespace precice
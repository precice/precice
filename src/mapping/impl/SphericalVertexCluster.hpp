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
   * matrices and computes the matrix decomposition direclty.
   *
   * @param[in] center Center of the particular partition
   * @param[in] radius radius of the partition
   * @param[in] parameter Shape parameter or support radius
   * @param[in] targetSize target size (cardinality) of the partition
   * @param[in] inputMesh the input mesh
   * @param[in] outputMesh the output mesh
   */
  SphericalVertexCluster(mesh::Vertex      center,
                         double            radius,
                         double            parameter,
                         std::vector<bool> deadAxis,
                         Polynomial        polynomial,
                         unsigned int      targetSize,
                         mesh::PtrMesh     inputMesh,
                         mesh::PtrMesh     outputMesh);

  /// Invalidates and erases the data structures the cluster holds
  void clear();

  /// Evaluates a conservative mapping and agglomerates the result in the given output data
  void mapConservative(mesh::PtrData inputData, mesh::PtrData outputData) const;

  /// Evaluates a consistent mapping and agglomerates the result in the given output data
  void mapConsistent(mesh::PtrData inputData, mesh::PtrData outputData) const;

  /// set the normalized weight in the normalizedWeight data structure
  void setNormalizedWeight(double normalizedWeight, VertexID vertexID);

  /// Compute the weight for a given vertex
  double computeWeight(const mesh::Vertex &v) const;

  /// Indicate whether a vertex lies within the partition or not
  bool isVertexInside(const mesh::Vertex &v) const;

  /// Number of input vertices this partition operates on
  unsigned int getNumberOfInputVertices() const;

  /// The center coordinate of this partition
  std::array<double, 3> getCenterCoords() const;

  /// Indicates whether there are valid output points for this partition
  bool empty() const;

private:
  precice::logging::Logger _log{"mapping::SphericalVertexCluster"};

  bool _hasComputedMapping = false;

  /// SphericalVertexCluster center
  mesh::Vertex _center;

  /// SphericalVertexCluster radius
  const double _radius;

  /// vector containing the normalized weights used for the output mesh data
  /// (consistent mapping) or input mesh data (conservative data)
  Eigen::VectorXd _normalizedWeights;

  // Stores the global IDs of the vertices so that we can apply a binary
  // search in order to query specific objects
  boost::container::flat_set<VertexID> _inputIDs;
  boost::container::flat_set<VertexID> _outputIDs;

  /// The RBF solver
  RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> _rbfSolver;

  /// Polynomial treatment in the RBF solver
  Polynomial _polynomial;

  //TODO: Can probably be removed
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// The weighting function
  CompactPolynomialC2 _weightingFunction;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::SphericalVertexCluster(
    mesh::Vertex      center,
    double            radius,
    double            parameter,
    std::vector<bool> deadAxis,
    Polynomial        polynomial,
    unsigned int      targetSize,
    mesh::PtrMesh     inputMesh,
    mesh::PtrMesh     outputMesh)
    : _center(center), _radius(radius), _polynomial(polynomial), _basisFunction(parameter), _weightingFunction(radius)
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
  _inputIDs.insert(inIDs.begin(), inIDs.end());
  _outputIDs.insert(outIDs.begin(), outIDs.end());

  // If the cluster is empty, we return immediately
  if (empty()) {
    return;
  }

  PRECICE_ASSERT(inIDs.size() > 0, "The source partition is empty whereas the target partition is non-empty.", inIDs.size(), outIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster input size: {}", inIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster output size: {}", outIDs.size());

  // Otherwise or system is undertermined
  // if (inIDs.size() < dimension + 1)
  //   _polynomial = Polynomial::OFF;

  // Construct the solver
  _rbfSolver          = RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{_basisFunction, *inputMesh.get(), _inputIDs, *outputMesh.get(), _outputIDs, deadAxis, _polynomial};
  _hasComputedMapping = true;
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

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::setNormalizedWeight(double normalizedWeight, VertexID id)
{
  PRECICE_ASSERT(_outputIDs.size() > 0);
  if (_normalizedWeights.size() == 0)
    _normalizedWeights.resize(_outputIDs.size());

  auto localID = _outputIDs.index_of(_outputIDs.find(id));
  PRECICE_ASSERT(_outputIDs.contains(id), id);
  PRECICE_ASSERT(localID < _normalizedWeights.size(), localID, _normalizedWeights.size());
  PRECICE_ASSERT(normalizedWeight > 0);

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
  // outputIDs and input as for conservative mappings in and output are swapped in terms of the mesh
  Eigen::VectorXd in(_outputIDs.size());
  in.setZero();

  // The result can directly be written into the global data structures
  Eigen::VectorXd result;
  // For every data component, perform mapping

  // First, extract the relevant input data from the global input data and store
  // it in a contiguous array
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Fill input from input data values (last polyparams entries remain zero)
    for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
      const auto dataIndex = *(_outputIDs.nth(i));
      in[i]                = localInData[dataIndex * nComponents + c] * _normalizedWeights[i];
    }
    result = _rbfSolver.solveConservative(in, _polynomial);

    PRECICE_ASSERT(result.size() == _inputIDs.size());

    // Now accumulate the result into our global output data
    for (unsigned int i = 0; i < _inputIDs.size(); ++i) {
      const auto dataIndex = *(_inputIDs.nth(i));
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

  Eigen::VectorXd in(_rbfSolver.getEvaluationMatrix().cols()); // rows == n
  in.setZero();

  // For every data component, perform mapping

  // First, extract the relevant input data from the global input data and store
  // it in a contiguous array
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Fill input from input data values (last polyparams entries remain zero)
    for (unsigned int i = 0; i < _inputIDs.size(); i++) {
      const auto dataIndex = *(_inputIDs.nth(i));
      in[i]                = localInData[dataIndex * nComponents + c];
    }

    auto result = _rbfSolver.solveConsistent(in, _polynomial);

    for (unsigned int i = 0; i < _outputIDs.size(); ++i) {
      const auto dataIndex = *(_outputIDs.nth(i));
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
bool SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::isVertexInside(const mesh::Vertex &v) const
{
  // TODO: We don't need to reduce the dead coordinates here as the values should reduce anyway
  auto res = computeSquaredDifference(_center.rawCoords(), v.rawCoords(), {{true, true, true}});
  return res < math::pow_int<2>(_radius);
}

} // namespace mapping
} // namespace precice

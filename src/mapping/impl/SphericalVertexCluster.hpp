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
 * @brief TODO
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class SphericalVertexCluster {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] dimension Dimensionality of the meshes
   * @param[in] center Center of the particular partition
   * @param[in] radius radius of the partition
   * @param[in] parameter Shape parameter or support radius
   * @param[in] targetSize target size (cardinality) of the partition
   * @param[in] inputMesh the input mesh
   * @param[in] outputMesh the output mesh
   */
  SphericalVertexCluster(int               dimension,
                         mesh::Vertex      center,
                         double            radius,
                         double            parameter,
                         std::vector<bool> deadAxis,
                         Polynomial        polynomial,
                         unsigned int      targetSize,
                         mesh::PtrMesh     inputMesh,
                         mesh::PtrMesh     outputMesh);

  /// Removes a computed mapping.
  void clear();

  /// Execute a consistent mapping
  void mapConsistent(mesh::PtrData inputData, mesh::PtrData outputData);

  /// Execute a consistent mapping
  void mapConservative(mesh::PtrData inputData, mesh::PtrData outputData);

  /// @brief set the normalized weight in the normalizedWeight data structure
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
  bool isEmpty() const;

private:
  precice::logging::Logger _log{"mapping::SphericalVertexCluster"};

  bool _hasComputedMapping = false;
  bool _emptyPartition     = false;

  /// SphericalVertexCluster center
  mesh::Vertex _center;

  /// SphericalVertexCluster radius
  const double _radius;

  /// vector containing the normalized weights (only for conservative mappings)
  Eigen::VectorXd normalizedWeights;

  // Stores the global IDs of the vertices so that we can apply a binary
  // search in order to query specific objects
  boost::container::flat_set<VertexID> inputIDs;
  boost::container::flat_set<VertexID> outputIDs;

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
    int               dimension,
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
  PRECICE_ASSERT(_polynomial != Polynomial::ON, "Integrated polynomial is not supported for PoU.")
  PRECICE_ASSERT(deadAxis.size() == dimension);
  PRECICE_DEBUG("Center coordinates: {}", _center.getCoords());
  PRECICE_DEBUG("SphericalVertexCluster radius: {}", _radius);

  // Get vertices to be mapped
  auto outIDs = outputMesh->index().getVerticesInsideBox(center, radius);
  // Constructing the partition when we don't have evaluation points is pointless
  auto inIDs = inputMesh->index().getVerticesInsideBox(center, radius);

  // If we have only very few vertices (10% of the target), we consider this partition as empty
  // TODO: The second empty criterion is potentially not suitable for conservative mappings
  _emptyPartition = outIDs.size() == 0 || inIDs.size() == 0;

  if (_emptyPartition) {
    return;
  }

  PRECICE_ASSERT(inIDs.size() > 0, "The source partition is empty whereas the target partition is non-empty.", inIDs.size(), outIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster input size: {}", inIDs.size());
  PRECICE_DEBUG("SphericalVertexCluster output size: {}", outIDs.size());

  // Otherwise or system is undertermined
  // if (inIDs.size() < dimension + 1)
  //   _polynomial = Polynomial::OFF;

  // Construct the solver
  inputIDs.insert(inIDs.begin(), inIDs.end());
  outputIDs.insert(outIDs.begin(), outIDs.end());
  _rbfSolver          = RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{_basisFunction, *inputMesh.get(), inputIDs, *outputMesh.get(), outputIDs, deadAxis, _polynomial};
  _hasComputedMapping = true;
}

template <typename RADIAL_BASIS_FUNCTION_T>
unsigned int SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::getNumberOfInputVertices() const
{
  return inputIDs.size();
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::array<double, 3> SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::getCenterCoords() const
{
  return _center.rawCoords();
}

template <typename RADIAL_BASIS_FUNCTION_T>
bool SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::isEmpty() const
{
  return _emptyPartition;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  inputIDs.clear();
  outputIDs.clear();
  _rbfSolver.clear();
  _hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::setNormalizedWeight(double normalizedWeight, VertexID id)
{
  PRECICE_ASSERT(outputIDs.size() > 0);
  if (normalizedWeights.size() == 0)
    normalizedWeights.resize(outputIDs.size());

  auto localID = outputIDs.index_of(outputIDs.find(id));
  PRECICE_ASSERT(outputIDs.contains(id), id);
  PRECICE_ASSERT(localID < normalizedWeights.size(), localID, normalizedWeights.size());
  PRECICE_ASSERT(normalizedWeight > 0);

  normalizedWeights[localID] = normalizedWeight;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConservative(mesh::PtrData inputData, mesh::PtrData outputData)
{
  // Empty partitions should not be stored at all
  PRECICE_ASSERT(!_emptyPartition);

  // Serial case
  PRECICE_ASSERT(_hasComputedMapping);
  const unsigned int nComponents = inputData->getDimensions();

  PRECICE_ASSERT(normalizedWeights.size() == outputIDs.size());
  const auto &localInData = inputData->values();

  // TODO: We can probably reduce the temporary allocations here
  // outputIDs and input as for conservative mappings in and output are swapped in terms of the mesh
  Eigen::VectorXd in(outputIDs.size());
  in.setZero();

  // The result can directly be written into the global data structures
  Eigen::VectorXd result;
  // For every data component, perform mapping

  // First, extract the relevant input data from the global input data and store
  // it in a contiguous array
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Fill input from input data values (last polyparams entries remain zero)
    for (unsigned int i = 0; i < outputIDs.size(); ++i) {
      const auto dataIndex = *(outputIDs.nth(i));
      in[i]                = localInData[dataIndex * nComponents + c] * normalizedWeights[i];
    }
    result = _rbfSolver.solveConservative(in, _polynomial);

    PRECICE_ASSERT(result.size() == inputIDs.size());

    // Now accumulate the result into our global output data
    for (unsigned int i = 0; i < inputIDs.size(); ++i) {
      const auto dataIndex = *(inputIDs.nth(i));
      outputData->values()[dataIndex * nComponents + c] += result(i);
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void SphericalVertexCluster<RADIAL_BASIS_FUNCTION_T>::mapConsistent(mesh::PtrData inputData, mesh::PtrData outputData)
{
  // Empty partitions should not be stored at all
  PRECICE_ASSERT(!_emptyPartition);

  // Serial case
  PRECICE_ASSERT(_hasComputedMapping);
  const unsigned int nComponents = inputData->getDimensions();

  const auto &localInData = inputData->values();

  Eigen::VectorXd in(_rbfSolver.getEvaluationMatrix().cols()); // rows == n
  in.setZero();

  // For every data component, perform mapping

  // First, extract the relevant input data from the global input data and store
  // it in a contiguous array
  for (unsigned int c = 0; c < nComponents; ++c) {
    // Fill input from input data values (last polyparams entries remain zero)
    for (unsigned int i = 0; i < inputIDs.size(); i++) {
      const auto dataIndex = *(inputIDs.nth(i));
      in[i]                = localInData[dataIndex * nComponents + c];
    }

    auto result = _rbfSolver.solveConsistent(in, _polynomial);

    for (unsigned int i = 0; i < outputIDs.size(); ++i) {
      const auto dataIndex = *(outputIDs.nth(i));
      outputData->values()[dataIndex * nComponents + c] += result(i) * normalizedWeights[i];
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

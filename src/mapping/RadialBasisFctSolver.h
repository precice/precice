#include "impl/BasisFunctions.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Event.hpp"

namespace precice {
extern bool syncMode;

struct RadialBasisFctSolver {
public:
  RadialBasisFctSolver() = default;

  /// Assembles the system matrices and computes the decomposition of the interpolation matrix
  template <typename RADIAL_BASIS_FUNCTION_T>
  void computeDecomposition(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis);

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(const Eigen::VectorXd &inputData) const;

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData) const;

  // Clear all stored matrices
  void clear();

  Eigen::MatrixXd _matrixA;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qr;
};

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, std::vector<bool> deadAxis);

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis);

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver::computeDecomposition(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis)
{
  // First, assemble the interpolation matrix
  _qr = buildMatrixCLU(basisFunction, inputMesh, deadAxis).colPivHouseholderQr();

  PRECICE_CHECK(_qr.isInvertible(),
                "The interpolation matrix of the RBF mapping from mesh {} to mesh {} is not invertable. "
                "This means that the mapping problem is not well-posed. "
                "Please check if your coupling meshes are correct. Maybe you need to fix axis-aligned mapping setups "
                "by marking perpendicular axes as dead?",
                inputMesh.getName(), outputMesh.getName());

  // Second, assemble evaluation matrix
  _matrixA = buildMatrixA(basisFunction, inputMesh, outputMesh, deadAxis);
}

Eigen::VectorXd RadialBasisFctSolver::solveConservative(const Eigen::VectorXd &inputData) const
{
  // TODO: Avoid temporary allocations
  PRECICE_ASSERT(inputData.size() == _matrixA.rows());
  Eigen::VectorXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.size() == _matrixA.cols());
  return _qr.solve(Au);
}

Eigen::VectorXd RadialBasisFctSolver::solveConsistent(const Eigen::VectorXd &inputData) const
{
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  Eigen::VectorXd p = _qr.solve(inputData);
  PRECICE_ASSERT(p.size() == _matrixA.cols());
  return _matrixA * p;
}

void RadialBasisFctSolver::clear()
{
  _matrixA = Eigen::MatrixXd();
  _qr      = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>();
}

// ------- Non-Member Functions ---------

/// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
static inline Eigen::VectorXd reduceVector(
    const Eigen::VectorXd &  fullVector,
    const std::vector<bool> &deadAxis)
{
  int deadDimensions = 0;
  int dimensions     = deadAxis.size();
  for (int d = 0; d < dimensions; d++) {
    if (deadAxis[d])
      deadDimensions += 1;
  }
  PRECICE_ASSERT(dimensions > deadDimensions, dimensions, deadDimensions);
  Eigen::VectorXd reducedVector(dimensions - deadDimensions);
  int             k = 0;
  for (int d = 0; d < dimensions; d++) {
    if (not deadAxis[d]) {
      reducedVector[k] = fullVector[d];
      k++;
    }
  }
  return reducedVector;
}

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, std::vector<bool> deadAxis)
{
  const unsigned int inputSize  = inputMesh.vertices().size();
  const unsigned int dimensions = inputMesh.getDimensions();

  const unsigned int deadDimensions = std::count(deadAxis.begin(), deadAxis.end(), true);
  const unsigned int polyparams     = 1 + dimensions - deadDimensions;
  const unsigned int n              = inputSize + polyparams; // Add linear polynom degrees

  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);

  Eigen::MatrixXd matrixCLU(n, n);
  matrixCLU.setZero();

  for (int i = 0; i < inputSize; ++i) {
    for (int j = i; j < inputSize; ++j) {
      const auto &u   = inputMesh.vertices()[i].getCoords();
      const auto &v   = inputMesh.vertices()[j].getCoords();
      matrixCLU(i, j) = basisFunction.evaluate(reduceVector((u - v), deadAxis).norm());
    }

    const auto reduced = reduceVector(inputMesh.vertices()[i].getCoords(), deadAxis);

    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      matrixCLU(i, inputSize + 1 + dim) = reduced[dim];
    }
    matrixCLU(i, inputSize) = 1.0;
  }

  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();

  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis)
{
  const unsigned int inputSize  = inputMesh.vertices().size();
  const unsigned int outputSize = outputMesh.vertices().size();
  const unsigned int dimensions = inputMesh.getDimensions();

  const unsigned int deadDimensions = std::count(deadAxis.begin(), deadAxis.end(), true);
  const unsigned int polyparams     = 1 + dimensions - deadDimensions;
  const unsigned int n              = inputSize + polyparams; // Add linear polynom degrees

  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);
  matrixA.setZero();

  // Fill _matrixA with values
  for (int i = 0; i < outputSize; ++i) {
    for (int j = 0; j < inputSize; ++j) {
      const auto &u = outputMesh.vertices()[i].getCoords();
      const auto &v = inputMesh.vertices()[j].getCoords();
      matrixA(i, j) = basisFunction.evaluate(reduceVector((u - v), deadAxis).norm());
    }

    const auto reduced = reduceVector(outputMesh.vertices()[i].getCoords(), deadAxis);

    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      matrixA(i, inputSize + 1 + dim) = reduced[dim];
    }
    matrixA(i, inputSize) = 1.0;
  }
  return matrixA;
}
} // namespace precice
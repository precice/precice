#include <numeric>
#include "mapping/impl/BasisFunctions.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Event.hpp"

namespace precice {
namespace mapping {

class RadialBasisFctSolver {
public:
  /// Default constructor
  RadialBasisFctSolver() = default;

  /// Assembles the system matrices and computes the decomposition of the interpolation matrix
  template <typename RADIAL_BASIS_FUNCTION_T>
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis);

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(const Eigen::VectorXd &inputData) const;

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData) const;

  // Clear all stored matrices
  void clear();

  // Access to the evaluation matrix
  const Eigen::MatrixXd &getEvaluationMatrix() const;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qr;

  Eigen::MatrixXd _matrixA;
};

// ------- Non-Member Functions ---------

/// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
inline double computeSquaredDifference(
    const std::array<double, 3> &u,
    std::array<double, 3>        v,
    const std::array<bool, 3> &  activeAxis)
{
  // Substract the values and multiply out dead dimensions
  for (unsigned int d = 0; d < v.size(); ++d) {
    v[d] = (u[d] - v[d]) * static_cast<int>(activeAxis[d]);
  }
  return std::accumulate(v.begin(), v.end(), static_cast<double>(0.), [](auto &res, auto &val) { return res + val * val; });
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, std::array<bool, 3> activeAxis)
{
  const auto         inputSize      = inputMesh.vertices().size();
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  // Treat the 2D case as 3D case with dead axis
  const unsigned int dimensions = 3;
  const unsigned int polyparams = 1 + dimensions - deadDimensions;
  const auto         n          = inputSize + polyparams; // Add linear polynom degrees

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);

  Eigen::MatrixXd matrixCLU(n, n);
  matrixCLU.setZero();

  for (auto i : boost::irange<Eigen::Index>(0, inputSize)) {
    for (auto j : boost::irange<Eigen::Index>(i, inputSize)) {
      const auto &u                 = inputMesh.vertices()[i].rawCoords();
      const auto &v                 = inputMesh.vertices()[j].rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixCLU(i, j)               = basisFunction.evaluate(std::sqrt(squaredDifference));
    }

    const auto &u = inputMesh.vertices()[i].rawCoords();

    unsigned int k = 0;
    for (unsigned int d = 0; d < dimensions; ++d) {
      if (activeAxis[d]) {
        matrixCLU(i, inputSize + 1 + k) = u[d];
        ++k;
      }
    }
    matrixCLU(i, inputSize) = 1.0;
  }

  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();

  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::array<bool, 3> activeAxis)
{
  const auto         inputSize      = inputMesh.vertices().size();
  const auto         outputSize     = outputMesh.vertices().size();
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  // Treat the 2D case as 3D case with dead axis
  const unsigned int dimensions = 3;
  const unsigned int polyparams = 1 + dimensions - deadDimensions;
  const auto         n          = inputSize + polyparams; // Add linear polynom degrees

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);
  matrixA.setZero();

  // Fill _matrixA with values

  for (auto i : boost::irange<Eigen::Index>(0, outputSize)) {
    for (auto j : boost::irange<Eigen::Index>(0, inputSize)) {
      const auto &u                 = outputMesh.vertices()[i].rawCoords();
      const auto &v                 = inputMesh.vertices()[j].rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixA(i, j)                 = basisFunction.evaluate(std::sqrt(squaredDifference));
    }

    const auto u = outputMesh.vertices()[i].rawCoords();

    unsigned int k = 0;
    for (unsigned int d = 0; d < dimensions; ++d) {
      if (activeAxis[d]) {
        matrixA(i, inputSize + 1 + k) = u[d];
        ++k;
      }
    }
    matrixA(i, inputSize) = 1.0;
  }
  return matrixA;
}

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctSolver::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis)
{
  // Convert dead axis vector into an active axis array so that we can handle the reduction more easily
  std::array<bool, 3> activeAxis({{false, false, false}});
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(), [](const auto ax) { return !ax; });
  // First, assemble the interpolation matrix
  _qr = buildMatrixCLU(basisFunction, inputMesh, activeAxis).colPivHouseholderQr();

  PRECICE_CHECK(_qr.isInvertible(),
                "The interpolation matrix of the RBF mapping from mesh {} to mesh {} is not invertable. "
                "This means that the mapping problem is not well-posed. "
                "Please check if your coupling meshes are correct. Maybe you need to fix axis-aligned mapping setups "
                "by marking perpendicular axes as dead?",
                inputMesh.getName(), outputMesh.getName());

  // Second, assemble evaluation matrix
  _matrixA = buildMatrixA(basisFunction, inputMesh, outputMesh, activeAxis);
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

const Eigen::MatrixXd &RadialBasisFctSolver::getEvaluationMatrix() const
{
  return _matrixA;
}
} // namespace mapping
} // namespace precice

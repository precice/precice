#include <numeric>
#include "mapping/config/MappingConfiguration.hpp"
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
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis, Polynomial polynomial);

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(Eigen::VectorXd &inputData, Polynomial polynomial) const;

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial) const;

  // Clear all stored matrices
  void clear();

  // Access to the evaluation matrix
  const Eigen::MatrixXd &getEvaluationMatrix() const;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qrMatrixC;

  // Decomposition of the polynomial
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qrMatrixQ;

  // TODO: Add max columns and rows in the template parameter
  Eigen::MatrixXd _matrixQ;
  Eigen::MatrixXd _matrixV;

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
  // @todo: this can be replaced by std::hypot when moving to C++17
  return std::accumulate(v.begin(), v.end(), static_cast<double>(0.), [](auto &res, auto &val) { return res + val * val; });
}

// Fill in the polynomial entries
inline void fillPolynomialEntries(Eigen::MatrixXd &matrix, const mesh::Mesh &mesh, Eigen::Index startIndex, std::array<bool, 3> activeAxis)
{
  // Loop over all vertices in the mesh
  for (auto i : boost::irange<Eigen::Index>(0, mesh.vertices().size())) {

    // 1. the constant contribution
    matrix(i, startIndex) = 1.0;

    // 2. the linear contribution
    const auto & u = mesh.vertices()[i].rawCoords();
    unsigned int k = 0;
    // Loop over all three space dimension and ignore dead axis
    for (unsigned int d = 0; d < activeAxis.size(); ++d) {
      if (activeAxis[d]) {
        PRECICE_ASSERT(matrix.rows() > i, matrix.rows(), i);
        PRECICE_ASSERT(matrix.cols() > startIndex + 1 + k, matrix.cols(), startIndex + 1 + k);
        matrix(i, startIndex + 1 + k) = u[d];
        ++k;
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  // Add linear polynom degrees if polynomial requires this
  const auto inputSize = inputMesh.vertices().size();
  const auto n         = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixCLU(n, n);
  matrixCLU.setZero();

  // Compute RBF matrix entries
  for (auto i : boost::irange<Eigen::Index>(0, inputSize)) {
    for (auto j : boost::irange<Eigen::Index>(i, inputSize)) {
      const auto &u                 = inputMesh.vertices()[i].rawCoords();
      const auto &v                 = inputMesh.vertices()[j].rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixCLU(i, j)               = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixCLU, inputMesh, inputSize, activeAxis);
  }
  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();
  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  const auto inputSize  = inputMesh.vertices().size();
  const auto outputSize = outputMesh.vertices().size();
  const auto n          = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);
  matrixA.setZero();

  // Compute RBF values for matrix A
  for (auto i : boost::irange<Eigen::Index>(0, outputSize)) {
    for (auto j : boost::irange<Eigen::Index>(0, inputSize)) {
      const auto &u                 = outputMesh.vertices()[i].rawCoords();
      const auto &v                 = inputMesh.vertices()[j].rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixA(i, j)                 = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixA, outputMesh, inputSize, activeAxis);
  }

  return matrixA;
}

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctSolver::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis, Polynomial polynomial)
{
  // Convert dead axis vector into an active axis array so that we can handle the reduction more easily
  std::array<bool, 3> activeAxis({{false, false, false}});
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(), [](const auto ax) { return !ax; });
  // First, assemble the interpolation matrix
  _qrMatrixC = buildMatrixCLU(basisFunction, inputMesh, activeAxis, polynomial).colPivHouseholderQr();

  PRECICE_CHECK(_qrMatrixC.isInvertible(),
                "The interpolation matrix of the RBF mapping from mesh {} to mesh {} is not invertable. "
                "This means that the mapping problem is not well-posed. "
                "Please check if your coupling meshes are correct. Maybe you need to fix axis-aligned mapping setups "
                "by marking perpendicular axes as dead?",
                inputMesh.getName(), outputMesh.getName());

  // Second, assemble evaluation matrix
  _matrixA = buildMatrixA(basisFunction, inputMesh, outputMesh, activeAxis, polynomial);

  // In case we deal with separated polynomials, we need dedicated matrices for the polynomial contribution
  if (polynomial == Polynomial::SEPARATE) {

    // 1. Allocate memory for these matrices
    // 4 = 1 + dimensions(3) = maximum number of polynomial parameters
    const unsigned int polyParams = 4 - std::count(activeAxis.begin(), activeAxis.end(), false);
    _matrixQ.resize(inputMesh.vertices().size(), polyParams);
    _matrixV.resize(outputMesh.vertices().size(), polyParams);

    // 2. fill the matrices: Q for the inputMesh, V for the outputMesh
    fillPolynomialEntries(_matrixQ, inputMesh, 0, activeAxis);
    fillPolynomialEntries(_matrixV, outputMesh, 0, activeAxis);

    // 3. compute decomposition
    _qrMatrixQ = _matrixQ.colPivHouseholderQr();
  }
}

Eigen::VectorXd RadialBasisFctSolver::solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial) const
{
  // TODO: Avoid temporary allocations
  // Au is equal to the eta in our PETSc implementation
  PRECICE_ASSERT(inputData.size() == _matrixA.rows());
  Eigen::VectorXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.size() == _matrixA.cols());

  // mu in the PETSc implementation
  Eigen::VectorXd out = _qrMatrixC.solve(Au);

  if (polynomial == Polynomial::SEPARATE) {
    Eigen::VectorXd epsilon = _matrixV.transpose() * inputData;
    PRECICE_ASSERT(epsilon.size() == _matrixV.cols());

    // epsilon = Q^T * mu - epsilon (tau in the PETSc impl)
    epsilon -= _matrixQ.transpose() * out;
    PRECICE_ASSERT(epsilon.size() == _matrixQ.cols());

    // out  = out - solveTranspose tau (sigma in the PETSc impl)
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
    out -= static_cast<Eigen::VectorXd>(_qrMatrixQ.transpose().solve(-epsilon));
#else
    // Backwards compatible version
    Eigen::VectorXd    sigma(_matrixQ.rows());
    const Eigen::Index nonzero_pivots = _qrMatrixQ.nonzeroPivots();

    if (nonzero_pivots == 0) {
      sigma.setZero();
    } else {
      Eigen::VectorXd c(_qrMatrixQ.colsPermutation().transpose() * (-epsilon));

      _qrMatrixQ.matrixQR().topLeftCorner(nonzero_pivots, nonzero_pivots).template triangularView<Eigen::Upper>().transpose().conjugate().solveInPlace(c.topRows(nonzero_pivots));

      sigma.topRows(nonzero_pivots) = c.topRows(nonzero_pivots);
      sigma.bottomRows(_qrMatrixQ.rows() - nonzero_pivots).setZero();

      sigma.applyOnTheLeft(_qrMatrixQ.householderQ().setLength(nonzero_pivots));
      out -= sigma;
    }
#endif
  }
  return out;
}

Eigen::VectorXd RadialBasisFctSolver::solveConsistent(Eigen::VectorXd &inputData, Polynomial polynomial) const
{
  Eigen::VectorXd res;
  // Solve polynomial QR and substract it form the input data
  if (polynomial == Polynomial::SEPARATE) {
    res = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * res);
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  Eigen::VectorXd p = _qrMatrixC.solve(inputData);
  PRECICE_ASSERT(p.size() == _matrixA.cols());
  Eigen::VectorXd out = _matrixA * p;

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * res);
  }

  return out;
}

void RadialBasisFctSolver::clear()
{
  _matrixA   = Eigen::MatrixXd();
  _qrMatrixC = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>();
}

const Eigen::MatrixXd &RadialBasisFctSolver::getEvaluationMatrix() const
{
  return _matrixA;
}
} // namespace mapping
} // namespace precice

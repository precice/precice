#pragma once

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
#include <numeric>
#include "mapping/MathHelper.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mesh/Mesh.hpp"
#include "profiling/Event.hpp"

namespace precice {
namespace mapping {

/// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
inline double computeSquaredDifference(
    const std::array<double, 3> &u,
    std::array<double, 3>        v,
    const std::array<bool, 3> &  activeAxis = {{true, true, true}})
{
  // Subtract the values and multiply out dead dimensions
  for (unsigned int d = 0; d < v.size(); ++d) {
    v[d] = (u[d] - v[d]) * static_cast<int>(activeAxis[d]);
  }
  // @todo: this can be replaced by std::hypot when moving to C++17
  return std::accumulate(v.begin(), v.end(), static_cast<double>(0.), [](auto &res, auto &val) { return res + val * val; });
}

/// given the active axis, computes sets the axis with the lowest spatial expansion to dead
template <typename IndexContainer>
constexpr void reduceActiveAxis(const mesh::Mesh &mesh, const IndexContainer &IDs, std::array<bool, 3> &axis)
{
  // make a pair of the axis and the difference
  std::array<std::pair<int, double>, 3> differences;

  // Compute the difference magnitude per direction
  for (std::size_t d = 0; d < axis.size(); ++d) {
    // Ignore dead axis here, i.e., apply the max value such that they are sorted on the last position(s)
    if (axis[d] == false) {
      differences[d] = std::make_pair<int, double>(d, std::numeric_limits<double>::max());
    } else {
      auto res = std::minmax_element(IDs.begin(), IDs.end(), [&](const auto &a, const auto &b) { return mesh.vertex(a).coord(d) < mesh.vertex(b).coord(d); });
      // Check if we are above or below the threshold
      differences[d] = std::make_pair<int, double>(d, std::abs(mesh.vertex(*res.second).coord(d) - mesh.vertex(*res.first).coord(d)));
    }
  }

  std::sort(differences.begin(), differences.end(), [](const auto &d1, const auto &d2) { return d1.second < d2.second; });
  // Disable the axis having the smallest expansion
  axis[differences[0].first] = false;
}

// Fill in the polynomial entries
template <typename IndexContainer>
inline void fillPolynomialEntries(Eigen::MatrixXd &matrix, const mesh::Mesh &mesh, const IndexContainer &IDs, Eigen::Index startIndex, std::array<bool, 3> activeAxis)
{
  // Loop over all vertices in the mesh
  for (const auto &i : IDs | boost::adaptors::indexed()) {

    // 1. the constant contribution
    matrix(i.index(), startIndex) = 1.0;

    // 2. the linear contribution
    const auto & u = mesh.vertex(i.value()).rawCoords();
    unsigned int k = 0;
    // Loop over all three space dimension and ignore dead axis
    for (unsigned int d = 0; d < activeAxis.size(); ++d) {
      if (activeAxis[d]) {
        PRECICE_ASSERT(matrix.rows() > i.index(), matrix.rows(), i.index());
        PRECICE_ASSERT(matrix.cols() > startIndex + 1 + k, matrix.cols(), startIndex + 1 + k);
        matrix(i.index(), startIndex + 1 + k) = u[d];
        ++k;
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                               std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  // Add linear polynom degrees if polynomial requires this
  const auto inputSize = inputIDs.size();
  const auto n         = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixCLU(n, n);

  // Required to fill the poly -> poly entries in the matrix, which remain otherwise untouched
  if (polynomial == Polynomial::ON) {
    matrixCLU.setZero();
  }

  // Compute RBF matrix entries
  auto         i_iter  = inputIDs.begin();
  Eigen::Index i_index = 0;
  for (; i_iter != inputIDs.end(); ++i_iter, ++i_index) {
    const auto &u       = inputMesh.vertex(*i_iter).rawCoords();
    auto        j_iter  = i_iter;
    auto        j_index = i_index;
    for (; j_iter != inputIDs.end(); ++j_iter, ++j_index) {
      const auto &v                 = inputMesh.vertex(*j_iter).rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixCLU(i_index, j_index)   = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixCLU, inputMesh, inputIDs, inputSize, activeAxis);
  }
  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();
  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                             const mesh::Mesh &outputMesh, const IndexContainer outputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  const auto inputSize  = inputIDs.size();
  const auto outputSize = outputIDs.size();
  const auto n          = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);

  // Compute RBF values for matrix A
  for (const auto &i : outputIDs | boost::adaptors::indexed()) {
    const auto &u = outputMesh.vertex(i.value()).rawCoords();
    for (const auto &j : inputIDs | boost::adaptors::indexed()) {
      const auto &v                 = inputMesh.vertex(j.value()).rawCoords();
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixA(i.index(), j.index()) = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixA, outputMesh, outputIDs, inputSize, activeAxis);
  }
  return matrixA;
}

// Variant operating on the Cholesky decopmosition
inline Eigen::VectorXd computeInverseDiagonal(Eigen::LLT<Eigen::MatrixXd> decMatrixC)
{
  // 1. Compute the diagonal entries of the inverse kernel matrix:
  // We already have the Cholesky decomposition. So instead of solving for the
  // kernel matrix directly, we invert the lower triangular matrix of the
  // decomposition:
  // using A^{-1} = (L^T)^{-1}L^{-1} enables the computation of the diagonal
  // entries of A^{-1} by evaluating the product above:
  // A^{-1}_{ii} = sum_{k=1}^n L^{-T}_{ik} L^{-1}_{ki}

  // 1a: Compute the inverse of the lower triangular matrix L
  // Eigen::MatrixXd L_inv = L.inverse(); is not supported by Eigen (linker errors)
  // However, Eigen provides triangular solver (LAPACK::trsm), which can be used
  // to solve L * Linv = I
  // Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(inSize, inSize);
  // _decMatrixC.matrixL().solveInPlace(L_inv);
  // which yields cubic complexity (BLAS level 3).

  // Here, we use our own implementation to compute the triangular inverse
  // (similar to LAPACK:trtri) more efficient, given that the RHS is also
  // triangular was unfortunately slower.

  // Solve L * Linv = I
  // Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(n, n);
  // decMatrixC.matrixL().solveInPlace(L_inv);
  Eigen::MatrixXd L_inv = utils::invertLowerTriangularBlockwise<double>(decMatrixC.matrixL());

  // 1b: Compute the diagonal elements of A^{-1} by evaluating (L^T)^{-1}L^{-1}
  Eigen::VectorXd inverseDiagonal = (L_inv.array().square().colwise().sum()).transpose();

  return inverseDiagonal;
}

// Variant operating on the QR decomposition
inline Eigen::VectorXd computeInverseDiagonal(Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decMatrixC)
{
  // 1. Compute the diagonal entries of the inverse kernel matrix:
  // We could use
  //     diag_inv_A= _decMatrixC.inverse().diagonal();
  // as Eigen offers this for the QR decomposition, but the inverse() call is
  // more expensive than it has to be (as we compute all entries of the inverse). On the
  // other hand, using non strictily-positive definite functions is less relevant anyway.
  // Still, let's try the following:

  // We already have the QR decomposition. So instead of solving for the
  // kernel matrix directly, make use of the following:
  // A^{-1} = R^{-1}Q^{-1}
  // Since Q is orthogonal, Q^{-1} = Q^T
  // R is upper triangular and we need to compute the inverse (using backwards substitution)

  // enables the computation of the diagonal
  // entries of A^{-1} by evaluating the product above:
  // A^{-1} = R^{-1} Q^T
  // A^{-1}_{ii} = sum_{k=1}^n R^{-1}_{ik} Q^{T}_{ki}

  Eigen::VectorXd    inverseDiagonal;
  const Eigen::Index n = decMatrixC.matrixR().cols();

  // 1a: Compute the inverse of the lower triangular matrix L
  // Solve R * Rinv = I
  Eigen::MatrixXd R_inv = Eigen::MatrixXd::Identity(n, n);
  decMatrixC.matrixR().template triangularView<Eigen::Upper>().solveInPlace(R_inv);
  Eigen::VectorXi P = decMatrixC.colsPermutation().indices();
  Eigen::MatrixXd Q = decMatrixC.householderQ();

  // Now evaluate, caution with the column permutation
  inverseDiagonal.resize(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    inverseDiagonal(P(i)) = R_inv.row(i) * Q.transpose().col(P(i));
  }
  return inverseDiagonal;
}

template<typename RBF_T>
inline Eigen::MatrixXd applyKernelToDistanceMatrix(const Eigen::MatrixXd &distanceMatrix, const Eigen::Index inputSize, const RBF_T &kernel)
{
  Eigen::MatrixXd matrixC = distanceMatrix;
  matrixC = matrixC.unaryExpr([&kernel] (double x) { return kernel.evaluate(x); }); // .block(0, 0, inputSize, inputSize)
  return matrixC;
}

inline double approximateConditionNumber(const Eigen::LLT<Eigen::MatrixXd> &choleskyDec)
{
  const Eigen::Index n = choleskyDec.matrixL().rows();

  double max_l = std::numeric_limits<double>::min();
  double min_l = std::numeric_limits<double>::max();

  for (Eigen::Index i = 0; i < n; i++) {
    const double lii = choleskyDec.matrixL()(i, i);
    min_l = std::min(lii, min_l);
    max_l = std::max(lii, min_l);
  }
  double condition = max_l * max_l / (min_l * min_l);
  if (condition < 0) condition = std::numeric_limits<double>::max();

  return condition;
}

inline double approximateConditionNumber(const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &qrDec)
{
  const Eigen::Index n = qrDec.matrixR().rows();

  double max_r = std::numeric_limits<double>::min();
  double min_r = std::numeric_limits<double>::max();

  for (Eigen::Index i = 0; i < n; i++) {
    const double rii = qrDec.matrixR()(i, i);
    min_r = std::min(rii, min_r);
    max_r = std::max(rii, min_r);
  }
  double condition = max_r / min_r; // TODO: correct?
  if (condition < 0) condition = std::numeric_limits<double>::max();

  return condition;
}

inline double computeRippaLOOCVerror(const Eigen::LLT<Eigen::MatrixXd> &decLLT, const Eigen::VectorXd &inputData)
{
  // Implementation of LOOCV according to Rippa(1999), DOI: 10.1023/a:1018975909870
  if (decLLT.info() != Eigen::ComputationInfo::Success) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const Eigen::Index n = inputData.size();
  const Eigen::VectorXd lambda = decLLT.solve(inputData);
  const double loocv = std::sqrt((lambda.array() / computeInverseDiagonal(decLLT).array()).array().square().sum() / n);

  return loocv;
}

}
}
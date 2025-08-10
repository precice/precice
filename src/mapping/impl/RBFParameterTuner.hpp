#pragma once

#include <functional>
#include <mapping/MathHelper.hpp>
#include <mapping/RadialBasisFctSolver.hpp>
#include <mapping/config/MappingConfigurationTypes.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>

namespace precice::mapping {

// Forward declaration of function found in <mapping/RadialBasisFctSolver.hpp>
// template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
// Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial);

struct Sample {
  double pos;
  double error;

  Sample(double pos, double error)
      : pos(pos), error(error)
  {
  }

  Sample()
      : pos(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::quiet_NaN())
  {
  }
};

template <typename Solver>
class RBFParameterTuner {

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

  using RBF_T = typename Solver::BASIS_FUNCTION_T; // TODO: better?

protected:

  Eigen::Index _inSize;
  double       _lowerBound;

  /*
  static constexpr bool rbfSupportsRadius()
  {
    return RBF_T::hasCompactSupport(); // TODO: necessary? better criterion?
  }

  static constexpr bool rbfUsesShapeParameter()
  {
    return std::is_same_v<RBF_T, Gaussian>; // TODO: necessary? better criterion?
  }
  */

  //using DecompositionType = std::conditional_t<RBF_T::isStrictlyPositiveDefinite(), Eigen::LLT<Eigen::MatrixXd>, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>;

public:
  virtual ~RBFParameterTuner() = default;

  RBFParameterTuner(const mesh::Mesh &inputMesh);
  virtual double optimize(const Solver &solver, const Eigen::VectorXd &inputData);

  //DecompositionType buildKernelDecomposition(double sampleRadius) const;

  //auto applyKernelToMatrix(const Eigen::Ref<const Eigen::MatrixXd> &matrixExpr, double sampleRadius) const;

  static double estimateMeshResolution(const mesh::Mesh &inputMesh);
};

template <typename Solver>
RBFParameterTuner<Solver>::RBFParameterTuner(const mesh::Mesh &inputMesh)
{
  _lowerBound = estimateMeshResolution(inputMesh);
  _inSize     = inputMesh.nVertices(); //TODO conversion?
}

template <typename Solver>
double RBFParameterTuner<Solver>::optimize(const Solver &solver, const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(false, "Not implemented.");
  return std::numeric_limits<double>::quiet_NaN();
}

/*
template <typename RBF_T>
const Eigen::MatrixXd &RBFParameterTuner<RBF_T>::getDistanceMatrix() const
{
  return _distanceMatrix;
}


template <typename Solver>
typename RBFParameterTuner<Solver>::DecompositionType RBFParameterTuner<Solver>::buildKernelDecomposition(double sampleRadius) const
{
  if constexpr (rbfSupportsRadius()) {
    double parameter = sampleRadius;
    if constexpr (rbfUsesShapeParameter()) {
      parameter = RBF_T::transformRadiusToShape(sampleRadius);
    }
    RBF_T kernel(parameter);
    // Check if kernel matrix was already initialized.
    if (_kernelMatrix.size() != _distanceMatrix.size()) {
      _kernelMatrix = _distanceMatrix;
    }
    // Apply kernel only to non-polynomial part of _distanceMatrix. The rest should remain unchanged.
    _kernelMatrix.block(0, 0, _inSize, _inSize) = _distanceMatrix.block(0, 0, _inSize, _inSize).unaryExpr([&kernel](double x) {
      return kernel.evaluate(x);
    });
    if constexpr (RBF_T::isStrictlyPositiveDefinite()) { // TODO: use inplace decomposition?
      return _kernelMatrix.llt();
    } else {
      return _kernelMatrix.colPivHouseholderQr();
    }
  }
  PRECICE_UNREACHABLE("RBF does not support a radius-initialization and was still used to instantiate an optimizer.");
}
*/


template <typename Solver>
double RBFParameterTuner<Solver>::estimateMeshResolution(const mesh::Mesh &inputMesh)
{
  constexpr int sampleSize = 3;

  const size_t       i0 = inputMesh.vertices().size() / 2;
  const mesh::Vertex x0 = inputMesh.vertices().at(i0);

  const std::vector<int> matches = inputMesh.index().getClosestVertices(x0.getCoords(), sampleSize);

  double h = 0;
  for (int i = 0; i < sampleSize; i++) {
    const mesh::Vertex xi = inputMesh.vertices().at(matches.at(i));
    h += std::sqrt(utils::computeSquaredDifference(xi.rawCoords(), x0.rawCoords()));
  }
  return h / sampleSize;
}

} // namespace precice::mapping

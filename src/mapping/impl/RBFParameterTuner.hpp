#pragma once

#include <functional>
#include <mapping/MathHelper.hpp>
#include <mapping/RadialBasisFctSolver.hpp>
#include <mapping/config/MappingConfigurationTypes.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>

namespace precice::mapping {

// Forward declaration of function found in <mapping/RadialBasisFctSolver.hpp>
template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial);

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

template <typename RBF_T>
class RBFParameterTuner {

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

protected:
  mutable Eigen::MatrixXd _kernelMatrix;
  Eigen::MatrixXd         _distanceMatrix;

  Eigen::Index _inSize;
  bool         _isInitialized;
  double       _lowerBound;

  static constexpr bool rbfSupportsRadius()
  {
    return RBF_T::hasCompactSupport(); // TODO: necessary? better criterion?
  }

  static constexpr bool rbfUsesShapeParameter()
  {
    return std::is_same_v<RBF_T, Gaussian>; // TODO: necessary? better criterion?
  }

  using DecompositionType = std::conditional_t<RBF_T::isStrictlyPositiveDefinite(), Eigen::LLT<Eigen::MatrixXd>, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>;

public:
  virtual ~RBFParameterTuner() = default;

  template <typename IndexContainer>
  RBFParameterTuner(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis);

  virtual double optimize(const Eigen::VectorXd &inputData);

  const Eigen::MatrixXd &getDistanceMatrix() const;
  DecompositionType      buildKernelDecomposition(double sampleRadius) const;

  using UnaryExprType = Eigen::CwiseUnaryOp<std::function<double(double)>, const Eigen::Ref<Eigen::MatrixXd>>;
  auto applyKernelToMatrix(const Eigen::Ref<const Eigen::MatrixXd> &matrixExpr, double sampleRadius) const;

  static double estimateMeshResolution(const mesh::Mesh &inputMesh);
  static double getMinBoundSize(const mesh::Mesh &inputMesh);
};

template <typename RBF_T>
template <typename IndexContainer>
RBFParameterTuner<RBF_T>::RBFParameterTuner(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis)
    : _kernelMatrix(Eigen::MatrixXd(0, 0))
{
  _lowerBound     = estimateMeshResolution(inputMesh);
  _inSize         = inputIDs.size();
  _distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
  _isInitialized  = true;
}

template <typename RBF_T>
double RBFParameterTuner<RBF_T>::optimize(const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(false, "Not implemented.");
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename RBF_T>
const Eigen::MatrixXd &RBFParameterTuner<RBF_T>::getDistanceMatrix() const
{
  return _distanceMatrix;
}

template <typename RBF_T>
typename RBFParameterTuner<RBF_T>::DecompositionType RBFParameterTuner<RBF_T>::buildKernelDecomposition(double sampleRadius) const
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

template <typename RBF_T>
auto RBFParameterTuner<RBF_T>::applyKernelToMatrix(const Eigen::Ref<const Eigen::MatrixXd> &matrixExpr, double sampleRadius) const
{
  std::function<double(double)> expr;

  if constexpr (rbfSupportsRadius()) {
    double parameter = sampleRadius;
    if constexpr (rbfUsesShapeParameter()) {
      parameter = RBF_T::transformRadiusToShape(sampleRadius);
    }
    RBF_T kernel(parameter);
    expr = [&kernel](double x) { return kernel.evaluate(x); };
  } else {
    PRECICE_ASSERT(false, "Selected RBF does not support a radius.");
    expr = [&](double x) { return x; };
  }
  return matrixExpr.unaryExpr(expr);
}

template <typename RBF_T>
double RBFParameterTuner<RBF_T>::estimateMeshResolution(const mesh::Mesh &inputMesh)
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

template <typename RBF_T>
double RBFParameterTuner<RBF_T>::getMinBoundSize(const mesh::Mesh &inputMesh)
{
  // inputMesh.computeBoundingBox();
  const mesh::BoundingBox &boundingBox = inputMesh.getBoundingBox();

  const int dimensions = inputMesh.getDimensions();
  double    minLength  = std::numeric_limits<double>::max();

  for (int axis = 0; axis < dimensions; axis++) {
    minLength = std::min(minLength, boundingBox.getEdgeLength(axis));
  }
  return minLength / 2;
};

} // namespace precice::mapping

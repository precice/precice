#pragma once

#include <limbo/limbo.hpp>
#include <mapping/RadialBasisFctSolver.hpp>
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/InitSampling.hpp"
#include "mapping/impl/OptimizationParameters.hpp"

namespace precice {
namespace mapping {

// Forward declaration of function found in <mapping/RadialBasisFctSolver.hpp>
// TODO: move to different file?
template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                               std::array<bool, 3> activeAxis, Polynomial polynomial);

struct Sample {
  double pos;
  double error;
};

template <typename RBF_T>
class RBFParameterTuner {

protected:
  Eigen::MatrixXd _distanceMatrix;
  Eigen::MatrixXd _kernelMatrix;
  Eigen::Index _inSize;
  bool _isInitialized;

  static constexpr bool rbfSupportsRadius()
  {
    return RBF_T::hasCompactSupport() || std::is_same_v<RBF_T, ThinPlateSplines>; // TODO: necessary? better criterion?
  }

  static constexpr bool rbfUsesShapeParameter()
  {
    return std::is_same_v<RBF_T, Gaussian>; // TODO: necessary? better criterion?
  }

public:
  virtual ~RBFParameterTuner() = default;

  RBFParameterTuner()
    : _distanceMatrix(Eigen::MatrixXd(0, 0)), _kernelMatrix(Eigen::MatrixXd(0, 0)), _inSize(0), _isInitialized(false)
  { }

  virtual double optimize(const Eigen::VectorXd &inputData)
  {
    return 0;
  }

  const Eigen::MatrixXd &getDistanceMatrix() const
  {
    return _distanceMatrix;
  }

  Eigen::LLT<Eigen::MatrixXd> buildKernelLLT(double sampleRadius)
  {
    if constexpr (rbfSupportsRadius() && RBF_T::isStrictlyPositiveDefinite()) {
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
      _kernelMatrix.block(0, 0, _inSize, _inSize) = _distanceMatrix.block(0, 0, _inSize, _inSize).unaryExpr([&kernel] (double x) {
        return kernel.evaluate(x);
      });
      return _kernelMatrix.llt();
    }
    return Eigen::LLT<Eigen::MatrixXd>();
  }

  Eigen::MatrixXd applyKernelToMatrix(const Eigen::MatrixXd &matrix, double sampleRadius)
  {
    if constexpr (rbfSupportsRadius()) {
      double parameter = sampleRadius;
      if constexpr (rbfUsesShapeParameter()) {
        parameter = RBF_T::transformRadiusToShape(sampleRadius);
      }
      RBF_T kernel(parameter);
      return matrix.unaryExpr([&kernel] (double x) { return kernel.evaluate(x); });
    } else {
      PRECICE_ASSERT(false, "Selected RBF does not support a radius.");
      return matrix.unaryExpr([&] (double x) { return x; });
    }
  }

protected:
  static double estimateMeshResolution(const mesh::Mesh &inputMesh)
  {
    constexpr int sampleSize = 5;

    const size_t           i0      = inputMesh.vertices().size() / 2;
    const mesh::Vertex     x0      = inputMesh.vertices().at(i0);
    const std::vector<int> matches = inputMesh.index().getClosestVertices(x0.getCoords(), sampleSize);

    double h = 0;
    for (int i = 0; i < sampleSize; i++) {
      const mesh::Vertex xi = inputMesh.vertices().at(matches.at(i));
      h += std::sqrt(utils::computeSquaredDifference(xi.rawCoords(), x0.rawCoords()));
    }
    return h / sampleSize;
  }
};

template<typename RBF_T>
class RBFParameterTunerSimple : public RBFParameterTuner<RBF_T> {

  double _lowerBound;
  std::vector<Sample> _samples;

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

public:
  RBFParameterTunerSimple()
      : _lowerBound(std::numeric_limits<double>::quiet_NaN())
  {
    PRECICE_ASSERT(RBF_T::isStrictlyPositiveDefinite(), "Non SPD RBFs are currently not supported by this optimizer");
    PRECICE_ASSERT(this->rbfSupportsRadius(), "RBF is not supported by this optimizer, as it does not accept a support-radius."
                                              "Currently supported: Compactly supported RBFs, Thin Plate Splines and Gaussians.");
  }

  template <typename IndexContainer>
  void initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool,3> &activeAxis)
  {
    _lowerBound   = this->estimateMeshResolution(inputMesh);
    this->_inSize = inputIDs.size();
    this->_distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
    this->_isInitialized  = true;
  }

  double golden_section_search(Sample lowerSample, Sample upperSample, const Eigen::VectorXd &inputData)
  {
    constexpr int    MAX_ITERATIONS = 10;
    constexpr double TOLERANCE = 1e-10;
    constexpr double PHI = (1.0 + std::sqrt(5)) / 2.0;

    Sample sample1;
    Sample sample2;

    sample1.pos = upperSample.pos + (lowerSample.pos - upperSample.pos) / PHI;
    sample2.pos = lowerSample.pos + (upperSample.pos - lowerSample.pos) / PHI;

    Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sample1.pos);
    sample1.error = utils::computeRippaLOOCVerror(llt, inputData); // quiet NaN if llt unsuccessful
    llt = this->buildKernelLLT(sample2.pos);
    sample2.error = utils::computeRippaLOOCVerror(llt, inputData);

    fmt::println("Start golden-section search:");

    for (int i = 0; i < MAX_ITERATIONS && std::abs(upperSample.error - lowerSample.error) > TOLERANCE; i++) {

      fmt::println("i={}: [ ({},{}) ({},{}) ({},{}) ({},{})], eps={}",
        i,
        lowerSample.pos, lowerSample.error,
        sample1.pos, sample1.error,
        sample2.pos, sample2.error,
        upperSample.pos, upperSample.error,
        std::abs(upperSample.error - lowerSample.error));

      if (sample1.error < sample2.error) {
        upperSample = sample2;
        sample2.pos = lowerSample.pos + (upperSample.pos - lowerSample.pos) / PHI;
        llt = this->buildKernelLLT(sample2.pos);
        sample2.error = utils::computeRippaLOOCVerror(llt, inputData);
      } else {
        lowerSample = sample1;
        sample1.pos = upperSample.pos + (lowerSample.pos - upperSample.pos) / PHI;
        llt = this->buildKernelLLT(sample1.pos);
        sample1.error = utils::computeRippaLOOCVerror(llt, inputData);
      }

      PRECICE_ASSERT(!std::isnan(sample1.error) && !std::isnan(sample2.error), "Optimization encountered NaN");

    }
    return (lowerSample.pos + upperSample.pos) / 2;
  }

  double optimize(const Eigen::VectorXd &inputData) override
  {
    PRECICE_ASSERT(this->_isInitialized);

    double increaseSize = 10;
    double sampleRadius = _lowerBound;
    double errorChange = 1e100;

    constexpr double MIN_INCREASE_SIZE      = 1.3;
    constexpr double MAX_NUMBER_OF_SAMPLES  = 10;
    constexpr double ERROR_CHANGE_TOLERANCE = 1e-10;

    // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
    while (errorChange > ERROR_CHANGE_TOLERANCE && increaseSize > MIN_INCREASE_SIZE && _samples.size() <= MAX_NUMBER_OF_SAMPLES) {

      // exponential increase of support-radius sampling position until failure
      while (errorChange > ERROR_CHANGE_TOLERANCE) {

        Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sampleRadius);
        if (llt.info() != Eigen::ComputationInfo::Success) {
          PRECICE_INFO("RBF tuner sample: rad={}, matrix decomposition failed", sampleRadius);
          PRECICE_CHECK(_samples.size() > 0, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
          sampleRadius = _samples.at(_samples.size() - 1).pos;
          break;
        }

        double rcond = utils::approximateReciprocalConditionNumber(llt);
        double error = utils::computeRippaLOOCVerror(llt, inputData);

        PRECICE_INFO("RBF tuner sample: rad={}, err={}, 1/cond={}", sampleRadius, error, rcond);

        if (rcond < 1e-13 || std::isnan(error)) {
          PRECICE_CHECK(_samples.size() > 0, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
          sampleRadius = _samples.at(_samples.size() - 1).pos;
          break;
        }
        errorChange = _samples.size() == 0 ? 1 : std::abs(_samples.at(_samples.size() - 1).error - error);
        PRECICE_INFO("RBF tuner sample: errorChange={}", errorChange); // TODO: not yet useful.

        _samples.push_back({sampleRadius, error});
        sampleRadius *= increaseSize;
      }
      increaseSize = std::sqrt(increaseSize);
      sampleRadius *= increaseSize;
    }

    size_t minIdx = 0;
    for (size_t i = 0; i < _samples.size(); i++) {
      if (_samples[i].error < _samples[minIdx].error) {
        minIdx = i;
      }
    }
    PRECICE_INFO("RBF optimization result: {}", _samples[minIdx].pos);
    return _samples[minIdx].pos;
  }
};


template<typename RBF_T>
class RBFParameterTunerBO : public RBFParameterTuner<RBF_T> {

public:
  using stop_t = boost::fusion::vector<limbo::stop::MaxIterations<OptimizationParameters>, limbo::stop::MaxPredictedValue<OptimizationParameters>>;
  using init_t = limbo::init::InitSampling<OptimizationParameters>;
  using acqui_t = limbo::acqui::EI<OptimizationParameters, limbo::model::GP<OptimizationParameters>>;
  using acqui_opt_t = limbo::opt::Cmaes<OptimizationParameters>;
  using ConfiguredBOOptimizer = limbo::bayes_opt::BOptimizer<OptimizationParameters, limbo::initfun<init_t>, limbo::acquifun<acqui_t>, limbo::acquiopt<acqui_opt_t>, limbo::stopcrit<stop_t>>;

  double _optimizedResult;
  double _lowerBound;
  int _n;

  // const mesh::Mesh &_inputMesh;
  Eigen::MatrixXd _distanceMatrix;

  // TODO: move


  // TODO: move
  double getMinBoundSize(const mesh::PtrMesh inputMesh) {
    inputMesh->computeBoundingBox(); // TODO: necessary?
    const mesh::BoundingBox &boundingBox = inputMesh->getBoundingBox();
    const int dimensions = inputMesh->getDimensions();
    double minLength = std::numeric_limits<double>::max();
    for (int axis = 0; axis < dimensions; axis++) {
      minLength = std::min(minLength, boundingBox.getEdgeLength(axis));
    }
    return minLength / 2;
  }

public:
  RBFParameterTunerBO()
      : _optimizedResult(-1), _lowerBound(-1), _n(-1), _distanceMatrix(Eigen::MatrixXd(0, 0))
  {
    // PRECICE_ASSERT(RBF_T::supportsRadiusInitialization());
    if (std::is_same_v<RBF_T, Gaussian>) { /* ... */
    }
  }



  template <typename IndexContainer>
  void initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool,3> &activeAxis)
  {
    PRECICE_ASSERT(RBF_T::isStrictlyPositiveDefinite());
    // PRECICE_ASSERT(isInitialized());
    _n = inputIDs.size();
    _lowerBound = this->estimateMeshResolution(inputMesh);
    _distanceMatrix = RadialBasisFctSolver<RBF_T>::buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
  }

  double optimize(const Eigen::VectorXd &inputData) {
    ConfiguredBOOptimizer _boptimizer;
    LOOCVEvaluator<RBF_T> _metricEvaluator(inputData, _distanceMatrix, _lowerBound);
    _boptimizer.optimize(_metricEvaluator);
    std::cout << "Best sample: " << _metricEvaluator.transformFromUnitToReal(_boptimizer.best_sample()(0)) << " - Best observation: " << _boptimizer.best_observation()(0) << std::endl;
    _optimizedResult = _metricEvaluator.transformFromUnitToReal(_boptimizer.best_sample()(0));
    return _optimizedResult;
  }

  double getOptimizedRadius() const
  {
    return _optimizedResult;
  }

  Eigen::MatrixXd getInterpolationMatrix()
  {
    RBF_T kernel(_optimizedResult);
    return applyKernelToDistanceMatrix(_distanceMatrix, _n, kernel);
  }



};

}
}
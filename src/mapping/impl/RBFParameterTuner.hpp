#pragma once

#include "mapping/MathHelper.hpp"
#include <mapping/RadialBasisFctSolver.hpp>
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/OptimizationParameters.hpp"
#include <limbo/limbo.hpp>
#include <libcmaes/cmaes.h>
#include <limbo/tools/macros.hpp>

namespace precice {
namespace mapping {

// Forward declaration of function found in <mapping/RadialBasisFctSolver.hpp>
// TODO: move to different file?
template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial);

struct Sample {
  double pos;
  double error;

  Sample(double pos, double error)
    : pos(pos), error(error)
  {}

  Sample()
    : pos(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::quiet_NaN())
  {}
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
    PRECICE_ASSERT(false, "Not implemented.");
    return std::numeric_limits<double>::quiet_NaN();
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
    double errorChange = 1e100; // TODO: unused

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

        _samples.emplace_back(sampleRadius, error);
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

  double optimize_bisection(const Eigen::VectorXd &inputData)
  {
    PRECICE_ASSERT(this->_isInitialized);

    double increaseSize = 10;
    double sampleRadius = _lowerBound;

    Sample lowerBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};
    Sample upperBound = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

    // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
    while (true) {
      sampleRadius *= increaseSize;

      Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sampleRadius);
      if (llt.info() != Eigen::ComputationInfo::Success) {
        PRECICE_INFO("RBF tuner sample: rad={}, matrix decomposition failed", sampleRadius);
        break;
      }

      double rcond = utils::approximateReciprocalConditionNumber(llt);
      double error = utils::computeRippaLOOCVerror(llt, inputData);

      PRECICE_INFO("RBF tuner sample: rad={}, err={}, 1/cond={}", sampleRadius, error, rcond);

      if (rcond < 1e-13 || std::isnan(error)) {
        PRECICE_CHECK(_samples.size() > 0, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
        break;
      }

      lowerBound = {sampleRadius, error};
    }
    upperBound = {sampleRadius, std::numeric_limits<double>::infinity()};


    constexpr double POS_TOLERANCE = 1.3;

    Sample centerSample;
    while (std::isnan(upperBound.error) && upperBound.pos < POS_TOLERANCE * lowerBound.pos) {
      centerSample.pos = (lowerBound.pos + upperBound.pos) / 2;

      Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sampleRadius);
      if (llt.info() != Eigen::ComputationInfo::Success) {
        centerSample.error = std::numeric_limits<double>::infinity();
      } else {
        centerSample.error = utils::computeRippaLOOCVerror(llt, inputData);
      }

      if (centerSample.error < upperBound.error || std::isinf(centerSample.error)) {
        upperBound = centerSample;
      }
    }

    return upperBound.pos;
  }
};


// ##################################################################################################################################################


namespace defaults {
struct init_initSampling {
  ///@ingroup init_defaults
  BO_PARAM(int, bins, 5);
};
} // namespace defaults


/** @ingroup init
  \rst
  Sampling which figures out the boundaries of out kernel.

  Parameter:
    - ``int bins`` (number of bins)
  \endrst
*/
template <typename Params>
struct InitSampling {

  mutable logging::Logger _log{"mapping::InitSampling"};

  template <typename StateFunction, typename AggregatorFunction, typename Opt> // hier: StateFunction = LOOCVEvaluator
  void operator()(const StateFunction &stateFn, const AggregatorFunction &, Opt &opt) const
  {
    std::vector<Sample> samples;

    double increaseSize = 10;
    double sampleRadius = stateFn.lower_bound;

    while (samples.size() < 2) {
      while (true) {

        const double error = stateFn(sampleRadius).value();
        const double rcond = stateFn.computeLastRCond();

        PRECICE_INFO("Bayes-Initializer: rad={}, err={}, 1/cond={}", sampleRadius, error, rcond);

        if (rcond < 1e-13 || std::isnan(error)) {
          break;
        }

        samples.emplace_back(sampleRadius, error);
        sampleRadius *= increaseSize;
      }
      increaseSize = std::sqrt(increaseSize);
      sampleRadius *= increaseSize;
    }

    stateFn.configureIntervalTransformation(samples);

    std::transform(samples.begin(), samples.end(), samples.begin(), [&stateFn](double sample) {
      return stateFn.transformSample(sample);
    });

    for (size_t i = 0; i < samples.size(); ++i) {
      opt.add_new_sample(limbo::tools::make_vector(samples[i].pos), limbo::tools::make_vector(samples[i].error));
    }
  }

};


template<typename RBF_T>
class RBFParameterTunerBO : public RBFParameterTuner<RBF_T> {

public:
  using stop_t = boost::fusion::vector<limbo::stop::MaxIterations<OptimizationParameters>, limbo::stop::MaxPredictedValue<OptimizationParameters>>;
  using init_t = InitSampling<OptimizationParameters>;
  using acqui_t = limbo::acqui::EI<OptimizationParameters, limbo::model::GP<OptimizationParameters>>;
  using acqui_opt_t = limbo::opt::Cmaes<OptimizationParameters>;
  using ConfiguredBayesOptimizer = limbo::bayes_opt::BOptimizer<OptimizationParameters, limbo::initfun<init_t>, limbo::acquifun<acqui_t>, limbo::acquiopt<acqui_opt_t>, limbo::stopcrit<stop_t>>;

  // number of input dimension (x.size())
  BO_PARAM(size_t, dim_in, 1);
  // number of dimensions of the result (res.size())
  BO_PARAM(size_t, dim_out, 1);

private:
  double _optimizedResult;
  double _lowerBound;
  double _upperBound;
  double _maxError;
  double _lastRCond;
  Eigen::VectorXd _inputData;

  Eigen::LLT<Eigen::MatrixXd> _llt;

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

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
      : _optimizedResult(-1), _lowerBound(std::numeric_limits<double>::infinity()),
        _upperBound(-std::numeric_limits<double>::infinity()),
        _maxError(-std::numeric_limits<double>::infinity()),
        _lastRCond(std::numeric_limits<double>::quiet_NaN())
  { }

  void configureIntervalTransformation(const std::vector<Sample> &samples)
  {
    for (const Sample &sample : samples) {
      _maxError   = std::max(_maxError, sample.error);
      _upperBound = std::min(_upperBound, sample.pos);
    }
  }

  bool isTransformationConfigured() const
  {
    return std::isfinite(_upperBound) && std::isfinite(_maxError);
  }

  double transformSamplePos(double pos) const
  {
    PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
    PRECICE_ASSERT(pos != 0, "Sample radius was 0 during transformation");
    return (std::log10(pos) - std::log10(_lowerBound)) / (std::log10(_upperBound) - std::log10(_lowerBound));
  }

  double backtransformSamplePos(double pos) const
  {
    PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
    return std::pow(10.0, pos * (std::log10(_upperBound) - std::log10(_lowerBound)) + std::log10(_lowerBound));
  }

  double transformSampleError(double error) const
  {
    PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
    return  error / _maxError;
  }

  double backtransformSampleError(double error) const
  {
    PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
    return error * _maxError;
  }

  Sample transformSample(Sample sample) const
  {
    sample.error = transformSampleError(sample.error);
    sample.pos = transformSamplePos(sample.pos);
    return sample;
  }

  template <typename IndexContainer>
  void initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool,3> &activeAxis)
  {
    PRECICE_ASSERT(RBF_T::isStrictlyPositiveDefinite());
    // PRECICE_ASSERT(isInitialized());
    _lowerBound = this->estimateMeshResolution(inputMesh);
    this->_inSize = inputIDs.size();
    this->_distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
  }

  double optimize(const Eigen::VectorXd &inputData) override
  {
    _inputData = inputData; // TODO: copy?
    ConfiguredBayesOptimizer bayesOptimizer;
    bayesOptimizer.optimize(*this);
    Sample bestSample = {backtransformSamplePos(bayesOptimizer.best_sample()(0)), backtransformSampleError(bayesOptimizer.best_observation()(0))};
    std::cout << "Best sample: " << bestSample.pos << " " << bestSample.error << std::endl;
    return bestSample.pos;
  }

  double getLastRCond() const
  {
    PRECICE_ASSERT(!std::isnan(_lastRCond), "No error has ben evaluated.");
    return _lastRCond;
  }

  Eigen::VectorXd operator()(const Eigen::VectorXd &x)
  {
    PRECICE_ASSERT(_inputData.size() > 0, "Input data not initialized.");

    double sampleRadius = isTransformationConfigured() ? transformSamplePos(x(0)) : x(0);

    _llt = this->buildKernelLLT(sampleRadius);
    double error = utils::computeRippaLOOCVerror(_llt, _inputData);
    error = isTransformationConfigured() ? transformSampleError(error) : error;
    _lastRCond = utils::approximateReciprocalConditionNumber(_llt);

    if (isTransformationConfigured()) {
      PRECICE_INFO("Bayes-Optimizer: rad={}, err={}", sampleRadius, error);
    }

    return limbo::tools::make_vector(error);
  }



};

}
}
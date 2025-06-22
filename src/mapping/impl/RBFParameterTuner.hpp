#pragma once

#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/OptimizationParameters.hpp"
#include "mapping/impl/InitSampling.hpp"
#include "mapping/impl/RBFMatrixOperations.hpp"
#include <limbo/limbo.hpp>

namespace precice {
namespace mapping {



class RBFParameterTuner {

public:
  virtual ~RBFParameterTuner() = default;
  struct Sample {
    double pos;
    double error;
  };

  virtual double optimize(const Eigen::VectorXd &inputData) { return 0; }

protected:
  double estimateMeshResolution(const mesh::Mesh &inputMesh) {
    constexpr int sampleSize = 5;
    const size_t i0 = inputMesh.vertices().size() / 2;
    const mesh::Vertex x0 = inputMesh.vertices().at(i0);
    const std::vector<int> matches = inputMesh.index().getClosestVertices(x0.getCoords(), sampleSize);
    double h = 0;
    for (int i = 0; i < sampleSize; i++) {
      const mesh::Vertex xi = inputMesh.vertices().at(matches.at(i));
      h += std::sqrt(computeSquaredDifference(xi.rawCoords(), x0.rawCoords()));
    }
    return h / sampleSize;
  }

};

template<typename RBF_T>
class RBFParameterTunerSimple : public RBFParameterTuner {


  double _lowerBound; // = this->estimateMeshResolution(inputMesh);
  double _upperBound;
  int _n;
  double _optimizedResult;
  Eigen::MatrixXd _distanceMatrix;
  std::vector<Sample> _samples;


  Eigen::LLT<Eigen::MatrixXd> computeLLT(double sampleRadius)
  {
    if constexpr (RBF_T::isStrictlyPositiveDefinite()) { // to make compiler happy (always true)
      RBF_T kernel(sampleRadius);
      return applyKernelToDistanceMatrix(_distanceMatrix, _n, kernel).llt();
    }
    return Eigen::LLT<Eigen::MatrixXd>();
  }

public:
  RBFParameterTunerSimple()
      : _lowerBound(0), _upperBound(std::numeric_limits<double>::infinity()), _n(-1), _optimizedResult(0), _distanceMatrix(Eigen::MatrixXd(0, 0))
  {
    PRECICE_ASSERT(RBF_T::isStrictlyPositiveDefinite()); //TODO: implementation for semidefinite kernels
    // PRECICE_ASSERT(RBF_T::supportsRadiusInitialization());
    if (std::is_same_v<RBF_T, Gaussian>) { /* ... */
    }
  }


  template <typename IndexContainer>
  void initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool,3> &activeAxis)
  {
    _n = inputIDs.size();
    _lowerBound = this->estimateMeshResolution(inputMesh);
    _distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
    fmt::println("init: _n={}", _n);
  }

  double optimize(const Eigen::VectorXd &inputData)
  {
    // PRECICE_ASSERT(isInitialized());

    // Find optimization interval

    double factor = 10;

    fmt::println("\nInitial sampling:\n");

    while (_samples.size() < 2) {
      double sampleRadius = _lowerBound;

      fmt::println("Start sampling with factor={} and rad={}", factor, sampleRadius);

      while (true) { // TODO: add max number of iterations

        Eigen::LLT<Eigen::MatrixXd> llt = computeLLT(sampleRadius);
        if (llt.info() != Eigen::ComputationInfo::Success) {
          fmt::println(" > Cholesky unsuccessful");
          break;
        }

        double condition = approximateConditionNumber(llt);        // evaluate cond at sampleRadius
        double error     = computeRippaLOOCVerror(llt, inputData); // evaluate error metric at sampleRadius

        fmt::println(" > rad={}, err={}, cond={}", sampleRadius, error, condition);

        if (condition > 1e8 || std::isnan(error) || std::isinf(error)) { // invalid sample
          break;
        }
        _samples.push_back({sampleRadius, error});
        sampleRadius *= factor;
      }
      factor = std::sqrt(factor);
    }
    _upperBound = _samples[_samples.size() - 1].pos;

    fmt::print("Resulting samples: [ ");
    for (auto &sample : _samples) {
      fmt::print("({},{}) ", sample.pos, sample.error);
    }
    fmt::println("]");

    // Reduce initial samples to a search interval

    size_t minIdx = 0;
    for (size_t i = 0; i < _samples.size() - 1; i++) {
      if (_samples[i].error < _samples[minIdx].error) {
        minIdx = i;
      }
    }

    Sample lowerSample;
    Sample upperSample;

    if (minIdx == 0) {
      lowerSample = _samples[minIdx];
      upperSample = _samples[minIdx + 1];
    } else if (minIdx == _samples.size() - 1) {
      lowerSample = _samples[minIdx - 1];
      upperSample = _samples[minIdx];
    } else {
      if (_samples[minIdx - 1].error < _samples[minIdx + 1].error) {
        lowerSample = _samples[minIdx - 1];
        upperSample = _samples[minIdx];
      } else {
        lowerSample = _samples[minIdx];
        upperSample = _samples[minIdx + 1];
      }
    }

    fmt::println("Search interval: [{}, {}]", lowerSample.pos, upperSample.pos);

    // Optimize using bisection
    // TODO: interval transformation?

    int maxOptimizationIterations = 10;
    double tolerance = 1e-10; // TODO: relative tolerance?
    constexpr double phi = (1.0 + std::sqrt(5)) / 2.0;

    Sample sample1;
    Sample sample2;

    sample1.pos = upperSample.pos + (lowerSample.pos - upperSample.pos) / phi;
    sample2.pos = lowerSample.pos + (upperSample.pos - lowerSample.pos) / phi;

    Eigen::LLT<Eigen::MatrixXd> llt = computeLLT(sample1.pos);
    sample1.error = computeRippaLOOCVerror(llt, inputData); // quiet NaN if llt unsuccessful
    llt = computeLLT(sample2.pos);
    sample2.error = computeRippaLOOCVerror(llt, inputData);

    fmt::println("Start golden-section search:");

    for (int i = 0; i < maxOptimizationIterations && std::abs(upperSample.error - lowerSample.error) < tolerance; i++) {

      fmt::println("i={}: [ ({},{}) ({},{}) ({},{}) ({},{})], i={}, eps={}",
        i,
        lowerSample.pos, lowerSample.error,
        sample1.pos, sample1.error,
        sample2.pos, sample2.error,
        upperSample.pos, upperSample.error,
        std::abs(upperSample.error - lowerSample.error));

      if (sample1.error < sample2.error) {
        upperSample = sample2;
        llt = computeLLT(sample2.pos);
        sample2.pos = lowerSample.pos + (upperSample.pos - lowerSample.pos) / phi;
        sample2.error = computeRippaLOOCVerror(llt, inputData);
      } else {
        lowerSample = sample1;
        llt = computeLLT(sample1.pos);
        sample1.pos = upperSample.pos + (lowerSample.pos - upperSample.pos) / phi;
        sample1.error = computeRippaLOOCVerror(llt, inputData);
      }

      PRECICE_ASSERT(!std::isnan(sample1.error) && !std::isnan(sample2.error), "Optimization encountered NaN"); // TODO: handle

    }

    _optimizedResult = (lowerSample.pos + upperSample.pos) / 2;
    fmt::println("Optimization result: {}", _optimizedResult);
    return _optimizedResult;
  }

  Eigen::MatrixXd getInterpolationMatrix()
  {
    RBF_T kernel(_optimizedResult);
    return applyKernelToDistanceMatrix(_distanceMatrix, _n, kernel);
  }
};


template<typename RBF_T>
class RBFParameterTunerBO : public RBFParameterTuner {

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
    _distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
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
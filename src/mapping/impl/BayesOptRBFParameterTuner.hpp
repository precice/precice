#pragma once

#include "mapping/impl/RBFParameterTuner.hpp"
#ifdef USE_LIBCMAES
  #include <libcmaes/cmaes.h>
#elif defined(USE_NLOPT)
  #include <nlopt.hpp>
#endif

#include <limbo/limbo.hpp>
#include <limbo/tools/macros.hpp>

namespace precice::mapping {

struct OptimizationParameters {

  struct bayes_opt_boptimizer : public limbo::defaults::bayes_opt_boptimizer {};

  // depending on which internal optimizer we use, we need to import different parameters
#ifdef USE_LIBCMAES
  struct opt_cmaes : public limbo::defaults::opt_cmaes {};
#elif defined(USE_NLOPT)
  struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {};
#else
  struct opt_gridsearch {
    BO_PARAM(int, bins, 20);
  };
#endif

  // enable / disable the writing of the result files
  struct bayes_opt_bobase {
    BO_PARAM(bool, stats_enabled, false);
    BO_PARAM(bool, bounded, true);
  };

  // no noise
  struct kernel : public limbo::defaults::kernel {
    BO_PARAM(double, noise, 1e-10);
  };

  struct kernel_maternfivehalves : public limbo::defaults::kernel_maternfivehalves {};

  // maximizing -(error / max error)
  // TODO: dependent on mesh width and kernel, requires knowledge of real interpolation error
  struct stop_maxpredictedvalue {
    BO_PARAM(double, ratio, -1e-15);
  };
  struct stop_maxiterations {
    BO_PARAM(int, iterations, 6);
  };

  // we use the default parameters for acqui_ucb
  // https://resibots.eu/limbo/api.html#acquisition-functions-acqui
  struct acqui_ucb {
    BO_PARAM(double, alpha, 1.5);
  };
  struct acqui_gpucb {
    BO_PARAM(double, delta, 0.2);
  };
  struct acqui_ei {
    BO_PARAM(double, jitter, 0.02); // ξ = jitter = 0.02
  };
};

struct InitSampling {

  mutable logging::Logger _log{"mapping::InitSampling"};

  template <typename StateFunction, typename AggregatorFunction, typename Opt> // here: StateFunction = RBFParameterTunerBO<RBF_T>
  void operator()(StateFunction &stateFn, const AggregatorFunction &, Opt &opt) const;
};

template <typename RBF_T>
class RBFParameterTunerBO : public RBFParameterTuner<RBF_T> {

  using DecompositionType = typename RBFParameterTuner<RBF_T>::DecompositionType;

public:
  using stop_t  = boost::fusion::vector<limbo::stop::MaxIterations<OptimizationParameters>, limbo::stop::MaxPredictedValue<OptimizationParameters>>;
  using init_t  = InitSampling;
  using acqui_t = limbo::acqui::EI<OptimizationParameters, limbo::model::GP<OptimizationParameters>>;

#ifdef USE_LIBCMAES
  using acquiopt_t = limbo::opt::Cmaes<OptimizationParameters>;
#elif defined(USE_NLOPT)
  using acquiopt_t = limbo::opt::NLOptNoGrad<OptimizationParameters, nlopt::GN_DIRECT_L_RAND>;
#else
  using acquiopt_t = limbo::opt::GridSearch<OptimizationParameters>;
#endif

  using ConfiguredBOptimizer = limbo::bayes_opt::BOptimizer<OptimizationParameters,
                                                            limbo::initfun<init_t>, limbo::acquifun<acqui_t>,
                                                            limbo::acquiopt<acquiopt_t>, limbo::stopcrit<stop_t>>;

  BO_PARAM(size_t, dim_in, 1);
  BO_PARAM(size_t, dim_out, 1);

private:


  // mutable variables because of const reference used by ConfiguredBOptimizer::optimize()
  mutable double _upperBound;
  mutable double _maxError;
  mutable double _lastRCond;

  Eigen::VectorXd _inputData;

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

public:
  template <typename IndexContainer>
  RBFParameterTunerBO(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis);

  double optimize(const Eigen::VectorXd &inputData) override;

  // Functions used by the limbo optimizer and InitSampling class

  Eigen::VectorXd operator()(const Eigen::VectorXd &x) const;

  double getLowerBound() const;
  double getLastRCond() const;

  void configureIntervalTransformation(const std::vector<Sample> &samples) const;
  bool isTransformationConfigured() const;

  double transformPos(double pos) const;
  double backtransformPos(double pos) const;
  double transformError(double error) const;
  double backtransformError(double error) const;
  Sample transformSample(Sample sample) const;
};

// Implementation:

template <typename StateFunction, typename AggregatorFunction, typename Opt>
void InitSampling::operator()(StateFunction &stateFn, const AggregatorFunction &, Opt &opt) const
{
  std::vector<Sample> samples;

  double increaseSize = 10;
  double sampleRadius = stateFn.getLowerBound();

  while (samples.size() < 2) {
    while (true) {

      const double error = stateFn(limbo::tools::make_vector(sampleRadius)).value();
      const double rcond = stateFn.getLastRCond();

      PRECICE_INFO("Bayes-Initializer: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

      samples.emplace_back(sampleRadius, error);

      if (rcond < 1e-13) {
        break;
      }
      sampleRadius *= increaseSize;
    }
    increaseSize = std::sqrt(increaseSize);
    sampleRadius *= increaseSize;
  }

  stateFn.configureIntervalTransformation(samples);

  std::transform(samples.begin(), samples.end(), samples.begin(), [&stateFn](Sample &sample) {
    return stateFn.transformSample(sample);
  });

  for (size_t i = 0; i < samples.size(); ++i) {
    opt.add_new_sample(limbo::tools::make_vector(samples[i].pos), limbo::tools::make_vector(samples[i].error));
  }
};

template <typename RBF_T> template <typename IndexContainer>
RBFParameterTunerBO<RBF_T>::RBFParameterTunerBO(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis)
    : RBFParameterTuner<RBF_T>(inputMesh, inputIDs, polynomial, activeAxis),
      _upperBound(-std::numeric_limits<double>::infinity()),
      _maxError(-std::numeric_limits<double>::infinity()),
      _lastRCond(std::numeric_limits<double>::quiet_NaN())
{
  PRECICE_ASSERT(this->rbfSupportsRadius(), "RBF is not supported by this optimizer, as it does not accept a support-radius."
                                            "Currently supported: Compactly supported RBFs and Gaussians.");
}

template <typename RBF_T>
void RBFParameterTunerBO<RBF_T>::configureIntervalTransformation(const std::vector<Sample> &samples) const
{
  for (const Sample &sample : samples) {
    _maxError   = std::max(_maxError, sample.error);
    _upperBound = std::max(_upperBound, sample.pos);
  }
}

template <typename RBF_T>
bool RBFParameterTunerBO<RBF_T>::isTransformationConfigured() const
{
  return std::isfinite(_upperBound) && std::isfinite(_maxError);
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::transformPos(double pos) const
{
  PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
  PRECICE_ASSERT(pos != 0, "Sample radius was 0 during transformation");
  return std::log10(pos / this->_lowerBound) / std::log10(_upperBound / this->_lowerBound);
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::backtransformPos(double pos) const
{
  PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
  return std::pow(10.0, pos * std::log10(_upperBound / this->_lowerBound) + std::log10(this->_lowerBound));
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::transformError(double error) const
{
  PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
  if (std::isnan(error)) {
    return -1;
  }
  return -error / _maxError;
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::backtransformError(double error) const
{
  PRECICE_ASSERT(isTransformationConfigured(), "Interval not configured.");
  return -error * _maxError;
}

template <typename RBF_T>
Sample RBFParameterTunerBO<RBF_T>::transformSample(Sample sample) const
{
  sample.error = transformError(sample.error);
  sample.pos   = transformPos(sample.pos);
  return sample;
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::getLowerBound() const
{
  PRECICE_ASSERT(this->_isInitialized);
  return this->_lowerBound;
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::optimize(const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(this->_isInitialized);
  _inputData = inputData; // TODO: avoid copy

  ConfiguredBOptimizer bayesOptimizer;
  bayesOptimizer.optimize(*this);

  double bestRadius = backtransformPos(bayesOptimizer.best_sample().value());
  double bestError  = backtransformError(bayesOptimizer.best_observation().value());
  PRECICE_INFO("Best sample: rad={:.4e}, err={:.4e}", bestRadius, bestError);

  return bestRadius;
}

template <typename RBF_T>
double RBFParameterTunerBO<RBF_T>::getLastRCond() const
{
  PRECICE_ASSERT(!std::isnan(_lastRCond), "No error has been evaluated.");
  return _lastRCond;
}

template <typename RBF_T>
Eigen::VectorXd RBFParameterTunerBO<RBF_T>::operator()(const Eigen::VectorXd &x) const
{
  PRECICE_ASSERT(_inputData.size() > 0, "Input data not initialized.");

  if (_maxError < 1e-15 && isTransformationConfigured()) {
    return limbo::tools::make_vector(0.0);
  }

  double sampleRadius = isTransformationConfigured() ? backtransformPos(x(0)) : x(0);

  const DecompositionType dec = this->buildKernelDecomposition(sampleRadius);
  if (!isTransformationConfigured()) {
    _lastRCond = utils::approximateReciprocalConditionNumber(dec);
  }

  double error = utils::computeRippaLOOCVerror(dec, _inputData);
  error        = isTransformationConfigured() ? transformError(error) : error;

  if (isTransformationConfigured()) {
    PRECICE_INFO("State function evaluation: f({:.4e})={:.4e}, [rad={:.4e}, err={:.4e}]", x(0), error, sampleRadius, backtransformError(error));
  }

  return limbo::tools::make_vector(error);
}

}

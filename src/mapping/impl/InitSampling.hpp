#pragma once
#include <Eigen/Core>

#include <limbo/tools/macros.hpp>

namespace limbo {
namespace defaults {
struct init_initSampling {
  ///@ingroup init_defaults
  BO_PARAM(int, bins, 5);
};
} // namespace defaults
namespace init {
/** @ingroup init
  \rst
  Sampling which figures out the boundaries of out kernel.

  Parameter:
    - ``int bins`` (number of bins)
  \endrst
*/
template <typename Params>
struct InitSampling {
  template <typename StateFunction, typename AggregatorFunction, typename Opt>
  void operator()(const StateFunction &seval, const AggregatorFunction &, Opt &opt) const
  {
    Eigen::VectorXd sample = Eigen::VectorXd::Constant(StateFunction::dim_in(), seval.lower_bound);

    std::vector<double> initSamples;
    std::vector<double> initObservations;

    double value         = seval(sample)(0); // Initial evaluation
    double currentSample = seval.lower_bound;

    // Add the initial sample if it is not NaN
    if (!std::isnan(value)) {
      initSamples.push_back(currentSample);
      initObservations.push_back(value);
    }

    // Increase sample by an order of magnitude
    while (true) {
      currentSample *= 10.0;
      sample = Eigen::VectorXd::Constant(StateFunction::dim_in(), currentSample);
      value  = seval(sample)(0);

      if (std::isnan(value)) {
        // Start bisection when a NaN is encountered
        double lower = initSamples.back(); // Last valid sample
        double upper = currentSample;

        for (int i = 0; i < 2; ++i) { // Perform at most 2 iterations
          currentSample = (lower + upper) / 2.0;
          sample        = Eigen::VectorXd::Constant(StateFunction::dim_in(), currentSample);
          value         = seval(sample)(0);

          if (!std::isnan(value)) {
            // Add valid sample and stop the bisection
            initSamples.push_back(currentSample);
            initObservations.push_back(value);
            break;
          } else {
            // Narrow down the interval
            upper = currentSample;
          }
        }
        break;
      } else {
        // Add valid sample
        initSamples.push_back(currentSample);
        initObservations.push_back(value);
      }
    }
    seval.setUpperLimit(initSamples.back());

    auto max = *std::max_element(initObservations.begin(), initObservations.end());
    auto min = *std::min_element(initObservations.begin(), initObservations.end());
    seval.setValues(min, max);

    std::transform(initObservations.begin(), initObservations.end(), initObservations.begin(), [&seval](double sample) {
      return seval.transformFromValues(sample);
    });

    std::transform(initSamples.begin(), initSamples.end(), initSamples.begin(), [&seval](double sample) {
      return seval.transformFromRealToUnit(sample);
    });
    for (int i = 0; i < initSamples.size(); ++i) {
      opt.add_new_sample(limbo::tools::make_vector(initSamples[i]), limbo::tools::make_vector(initObservations[i]));
    }
  }
};
} // namespace init
} // namespace limbo

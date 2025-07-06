#pragma once

#include <limbo/bayes_opt/boptimizer.hpp>
#include <limbo/init/grid_sampling.hpp>
#include <limbo/limbo.hpp>
#include <limbo/stop/max_predicted_value.hpp>

namespace precice {
struct OptimizationParameters {
  struct bayes_opt_boptimizer : public limbo::defaults::bayes_opt_boptimizer {
  };

// depending on which internal optimizer we use, we need to import different parameters
#ifdef USE_NLOPT
  struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {
  };
#elif defined(USE_LIBCMAES)
  struct opt_cmaes : public limbo::defaults::opt_cmaes {
  };
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

  struct kernel_maternfivehalves : public limbo::defaults::kernel_maternfivehalves {
  };

  // we use 10 random samples to initialize the algorithm
  struct init_initsampling {
    BO_PARAM(int, bins, 6);
  };
  struct stop_maxpredictedvalue {
    BO_PARAM(double, ratio, 0.9);
  };

  struct stop_maxiterations {
    BO_PARAM(int, iterations, 10);
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

} // namespace precice

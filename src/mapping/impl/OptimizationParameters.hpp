#pragma once
#include <limbo/bayes_opt/boptimizer.hpp>
#include <limbo/init/grid_sampling.hpp>
#include <limbo/limbo.hpp>
#include <limbo/stop/max_predicted_value.hpp>
#include <mapping/impl/InitSampling.hpp>
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

template <typename RADIAL_BASIS_FUNCTION_T>
struct LOOCVEvaluator {

  LOOCVEvaluator(const Eigen::VectorXd &in, const Eigen::MatrixXd &distance, double r)
      : lower_bound(r * 0.5), invalues(in), distanceMatrix(distance) {}

  // number of input dimension (x.size())
  BO_PARAM(size_t, dim_in, 1);
  // number of dimensions of the result (res.size())
  BO_PARAM(size_t, dim_out, 1);

  const double          lower_bound;
  mutable double        upper_bound = -1;
  const Eigen::VectorXd invalues;
  const Eigen::MatrixXd distanceMatrix;
  mutable double        lower_val = 1;
  mutable double        upper_val = 1;

  void setValues(double lower, double upper) const
  {
    lower_val = lower;
    upper_val = upper;
  }

  void setUpperLimit(double limit) const
  {
    upper_bound = limit;
  }

  double transformFromValues(double val) const
  {
    return val / std::abs(upper_val);
    /*
    if ((lower_val == upper_val) && (lower_val != 0)) {
      return val / std::abs(upper_val);
    } else {
      std::cout << "transform: " << (val - lower_val) / (upper_val - lower_val) << std::endl;
      return (val - lower_val) / (upper_val - lower_val);
    }
    */
  }

  double transformFromUnitToReal(double in) const
  {
    PRECICE_ASSERT(upper_bound > 0);
    return lower_bound + (upper_bound - lower_bound) * in;
  }

  double transformFromRealToUnit(double real) const
  {
    PRECICE_ASSERT(upper_bound > lower_bound);
    return (real - lower_bound) / (upper_bound - lower_bound);
  }

  // the function to be optimized
  Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
  {
    using RES_T = double;

    const Eigen::Index n = invalues.size();
    double x_eval = upper_bound == -1 ? x(0) : transformFromUnitToReal(x(0));

    std::unique_ptr<RADIAL_BASIS_FUNCTION_T> kernel;
    if constexpr (RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite()) {
      kernel = std::make_unique<RADIAL_BASIS_FUNCTION_T>(x_eval);
    } else {
      PRECICE_ASSERT(false, "Not supported.");
    }

    Eigen::LLT<Eigen::Matrix<RES_T, -1, -1>> dec = distanceMatrix.unaryExpr([&kernel](double x) { return static_cast<RES_T>(kernel->evaluate(x)); }).llt();

    // TODO: find a better value for failure
    if (dec.info() != Eigen::ComputationInfo::Success) {
      std::cout << std::fixed << std::setprecision(10) << "Parameter: " << x_eval << "  diverged." << std::endl;
      return limbo::tools::make_vector(std::numeric_limits<double>::quiet_NaN());
    }

    Eigen::Matrix<RES_T, -1, -1> L_inv  = utils::invertLowerTriangularBlockwise<RES_T>(dec.matrixL());
    Eigen::Vector<RES_T, -1>     incopy = invalues.unaryExpr([](double x) { return static_cast<RES_T>(x); });
    Eigen::Vector<RES_T, -1>     lambda = dec.solve(incopy);
    double                       loocv  = std::sqrt(static_cast<double>((lambda.array() / (L_inv.array().square().colwise().sum()).transpose().array()).array().square().sum()) / n);
    // Eigen::Vector<RES_T, -1> inv_diag = L_inv.array().square().colwise().sum().inverse();
    // double                   loocv = (((-0.5 * (lambda.array().square().array().colwise() * inv_diag.array())).array().colwise() - 0.5 * inv_diag.array().log().array()) - 0.5 * std::log(2 * M_PI)).colwise().sum().sum();
    std::cout << std::fixed << std::setprecision(10) << "Parameter: " << x_eval << "  with error: " << transformFromValues(-loocv) << std::endl;
    return limbo::tools::make_vector(transformFromValues(-loocv));
  }
};

} // namespace precice

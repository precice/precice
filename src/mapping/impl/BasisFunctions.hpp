#pragma once

#if defined(__NVCC__)

#include <cuda_runtime.h>
#define BOOST_PP_VARIADICS 1
#define PRECICE_HOST_DEVICE __host__ __device__
#define PRECICE_MEMORY_SPACE __device__
#define FMA fma

#elif defined(__HIPCC__)

// #define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define PRECICE_HOST_DEVICE __host__ __device__
#define PRECICE_MEMORY_SPACE __device__
#define FMA fma

#else

#define PRECICE_HOST_DEVICE
#define PRECICE_MEMORY_SPACE
#define FMA std::fma

#endif

PRECICE_MEMORY_SPACE const double NUMERICAL_ZERO_DIFFERENCE = 1.0e-14;

#include "logging/Logger.hpp"
#include "math/math.hpp"

namespace precice {
namespace mapping {

/**
 * @brief Wrapper struct that is used to transfer RBF-specific parameters to the GPU.
 *
 * Its parameters 1, 2 and 3 are filled individually, depending on the requirements
 * of each RBF (e.g., support radius, shape parameter, etc.). 
 * This parameter struct is handed to the GPU device kernel 
 * when assemblying the interpolation matrices
 * and its parameters used in every RBF evaluation.
 * Since there are at most 3 parameters for a RBF in the
 * current implementation, 3 parameters are required in this
 * struct, but ignored if the specific RBF does not require them all.
 *
 */
struct RadialBasisParameters {
  double parameter1;
  double parameter2;
  double parameter3;
};

/// Base class for RBF with compact support
struct CompactSupportBase {
  static constexpr bool hasCompactSupport()
  {
    return true;
  }
};

/// Base class for RBF without compact support
struct NoCompactSupportBase {
  static constexpr bool hasCompactSupport()
  {
    return false;
  }

  static constexpr double getSupportRadius()
  {
    return std::numeric_limits<double>::max();
  }
};

/**
 * @brief Base class for RBF functions to distinguish positive definite functions
 *
 * @tparam isDefinite
 */
template <bool isDefinite>
struct DefiniteFunction {
  static constexpr bool isStrictlyPositiveDefinite()
  {
    return isDefinite;
  }
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius^2 * log(radius).
 */
class ThinPlateSplines : public NoCompactSupportBase,
                         public DefiniteFunction<false> {
public:
  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    // We don't need to read any values from params since there is no need here
    return std::log(std::max(radius, NUMERICAL_ZERO_DIFFERENCE)) * math::pow_int<2>(radius);
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  RadialBasisParameters _params;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: sqrt(shape^2 + radius^2).
 */
class Multiquadrics : public NoCompactSupportBase,
                      public DefiniteFunction<false> {
public:
  explicit Multiquadrics(double c)
      : _cPow2(math::pow_int<2>(c))
  {
    _params.parameter1 = _cPow2;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double cPow2 = params.parameter1;
    return std::sqrt(cPow2 + math::pow_int<2>(radius));
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _cPow2;
  RadialBasisParameters _params;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: 1 / (shape^2 + radius^2).
 */
class InverseMultiquadrics : public NoCompactSupportBase,
                             public DefiniteFunction<true> {
public:
  explicit InverseMultiquadrics(double c)
      : _cPow2(math::pow_int<2>(c))
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::InverseMultiQuadrics"};
    PRECICE_CHECK(math::greater(c, 0.0),
                  "Shape parameter for radial-basis-function inverse multiquadric has to be larger than zero. Please update the \"shape-parameter\" attribute.");
#endif
    _params.parameter1 = _cPow2;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double cPow2 = params.parameter1;
    return 1.0 / std::sqrt(cPow2 + math::pow_int<2>(radius));
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double const _cPow2;

  RadialBasisParameters _params;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius.
 */
class VolumeSplines : public NoCompactSupportBase,
                      public DefiniteFunction<false> {
public:
  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    return std::abs(radius);
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  RadialBasisParameters _params;
};

/**
 * @brief Radial basis function with global and compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: exp(-1 * (shape * radius)^2).
 */
class Gaussian : public CompactSupportBase,
                 public DefiniteFunction<true> {
public:
  Gaussian(const double shape, const double supportRadius = std::numeric_limits<double>::infinity())
      : _shape(shape),
        _supportRadius(supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::Gaussian"};
    PRECICE_CHECK(math::greater(_shape, 0.0),
                  "Shape parameter for radial-basis-function gaussian has to be larger than zero. Please update the \"shape-parameter\" attribute.");
    PRECICE_CHECK(math::greater(_supportRadius, 0.0),
                  "Support radius for radial-basis-function gaussian has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif

    if (supportRadius < std::numeric_limits<double>::infinity()) {
      _deltaY = evaluate(supportRadius);
    }
    double threshold   = std::sqrt(-std::log(cutoffThreshold)) / shape;
    _supportRadius     = std::min(supportRadius, threshold);
    _params.parameter1 = _shape;
    _params.parameter2 = _supportRadius;
    _params.parameter3 = _deltaY;
  }

  double getSupportRadius() const
  {
    return _supportRadius;
  }

  double evaluate(const double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double shape         = params.parameter1;
    double supportRadius = params.parameter2;
    double deltaY        = params.parameter3;

    if (radius > supportRadius)
      return 0.0;
    return std::exp(-math::pow_int<2>(shape * radius)) - deltaY;
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  };

public:
  /// Below that value the function is supposed to be zero. Defines the support radius if not explicitly given
  static constexpr double cutoffThreshold = 1e-9;

private:
  double const _shape;

  /// Either explicitly set (from cutoffThreshold) or computed supportRadius
  double _supportRadius;

  double _deltaY = 0;

  RadialBasisParameters _params;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: 1 - 30*rn^2 - 10*rn^3 + 45*rn^4 - 6*rn^5 - 60*rn^3 * log(rn),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 * To work around the issue of log(0), the equation is formulated differently in the last term.
 */
class CompactThinPlateSplinesC2 : public CompactSupportBase,
                                  public DefiniteFunction<true> {
public:
  explicit CompactThinPlateSplinesC2(double supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::CompactThinPlateSplinesC2"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact thin-plate-splines c2 has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif
    _r_inv             = 1. / supportRadius;
    _params.parameter1 = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double       r_inv = params.parameter1;
    const double p     = radius * r_inv;
    if (p >= 1)
      return 0.0;
    return 1.0 - 30.0 * math::pow_int<2>(p) - 10.0 * math::pow_int<3>(p) + 45.0 * math::pow_int<4>(p) - 6.0 * math::pow_int<5>(p) - math::pow_int<3>(p) * 60.0 * std::log(std::max(p, NUMERICAL_ZERO_DIFFERENCE));
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^2,
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC0 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC0(double supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::CompactPolynomialC0"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c0 has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif
    _r_inv             = 1. / supportRadius;
    _params.parameter1 = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double       r_inv = params.parameter1;
    const double p     = radius * r_inv;
    if (p >= 1)
      return 0.0;
    return math::pow_int<2>(1.0 - p);
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^4 * ( 4rn + 1),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC2 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC2(double supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::CompactPolynomialC2"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c2 has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif

    _r_inv             = 1. / supportRadius;
    _params.parameter1 = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double       r_inv = params.parameter1;
    const double p     = radius * r_inv;
    if (p >= 1)
      return 0.0;
    return math::pow_int<4>(1.0 - p) * FMA(4, p, 1);
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^6 * ( 35 * (rn)^2 + 18rn + 3),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC4 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC4(double supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::CompactPolynomialC4"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c4 has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif

    _r_inv             = 1. / supportRadius;
    _params.parameter1 = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double       r_inv = params.parameter1;
    const double p     = radius * r_inv;
    if (p >= 1)
      return 0.0;
    return math::pow_int<6>(1.0 - p) * (35 * math::pow_int<2>(p) + FMA(18, p, 3));
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^8 * (32*rn^3 + 25*rn^2 + 8*rn + 1),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC6 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC6(double supportRadius)
  {
#if !defined(__NVCC__) || !defined(__HIPCC__)
    logging::Logger _log{"mapping::CompactPolynomialC6"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute.");
#endif
    _r_inv             = 1. / supportRadius;
    _params.parameter1 = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return operator()(radius, _params);
  }

  PRECICE_HOST_DEVICE inline double operator()(const double radius, const RadialBasisParameters params) const
  {
    double       r_inv = params.parameter1;
    const double p     = radius * r_inv;
    if (p >= 1)
      return 0.0;
    return math::pow_int<8>(1.0 - p) * (32.0 * math::pow_int<3>(p) + 25.0 * math::pow_int<2>(p) + FMA(8.0, p, 1.0));
  };

  const RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};
} // namespace mapping
} // namespace precice

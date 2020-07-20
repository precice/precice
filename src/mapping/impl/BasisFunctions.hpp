#pragma once

#include "logging/Logger.hpp"
#include "math/math.hpp"

namespace precice {
namespace mapping {

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
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius^2 * log(radius).
 */
class ThinPlateSplines : public NoCompactSupportBase {
public:
  double evaluate(double radius) const
  {
    double result = 0.0;
    if (math::greater(radius, 0.0)) {
      result = std::log(radius) * std::pow(radius, 2);
    }
    return result;
  }
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: sqrt(shape^2 + radius^2).
 */
class Multiquadrics : public NoCompactSupportBase {
public:
  explicit Multiquadrics(double c)
      : _cPow2(std::pow(c, 2)) {}

  double evaluate(double radius) const
  {
    return std::sqrt(_cPow2 + std::pow(radius, 2));
  }

private:
  double _cPow2;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: 1 / (shape^2 + radius^2).
 */
class InverseMultiquadrics : public NoCompactSupportBase {
public:
  explicit InverseMultiquadrics(double c)
      : _cPow2(std::pow(c, 2))
  {
    PRECICE_CHECK(math::greater(c, 0.0),
                  "Shape parameter for radial-basis-function inverse multiquadric has to be larger than zero. Please update the \"shape-parameter\" attribute.");
  }

  double evaluate(double radius) const
  {
    return 1.0 / std::sqrt(_cPow2 + std::pow(radius, 2));
  }

private:
  logging::Logger _log{"mapping::InverseMultiQuadrics"};

  double const _cPow2;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius.
 */
class VolumeSplines : public NoCompactSupportBase {
public:
  double evaluate(double radius) const
  {
    return std::abs(radius);
  }
};

/**
 * @brief Radial basis function with global and compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: exp(-1 * (shape * radius)^2).
 */
class Gaussian : public CompactSupportBase {
public:
  Gaussian(const double shape, const double supportRadius = std::numeric_limits<double>::infinity())
      : _shape(shape),
        _supportRadius(supportRadius)
  {
    PRECICE_CHECK(math::greater(_shape, 0.0),
                  "Shape parameter for radial-basis-function gaussian has to be larger than zero. Please update the \"shape-parameter\" attribute.");
    PRECICE_CHECK(math::greater(_supportRadius, 0.0),
                  "Support radius for radial-basis-function gaussian has to be larger than zero. Please update the \"support-radius\" attribute.");

    if (supportRadius < std::numeric_limits<double>::infinity()) {
      _deltaY = evaluate(supportRadius);
    }
    double threshold = std::sqrt(-std::log(cutoffThreshold)) / shape;
    _supportRadius   = std::min(supportRadius, threshold);
  }

  double getSupportRadius() const
  {
    return _supportRadius;
  }

  double evaluate(const double radius) const
  {
    if (radius > _supportRadius)
      return 0.0;
    else
      return std::exp(-std::pow(_shape * radius, 2.0)) - _deltaY;
  }

private:
  logging::Logger _log{"mapping::Gaussian"};

  /// Below that value the function is supposed to be zero. Defines the support radius if not explicitely given
  double const cutoffThreshold = 1e-9;

  double const _shape;

  /// Either explicitely set (from cutoffThreshold) or computed supportRadius
  double _supportRadius;

  double _deltaY = 0;
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
class CompactThinPlateSplinesC2 : public CompactSupportBase {
public:
  explicit CompactThinPlateSplinesC2(double supportRadius)
      : _r(supportRadius)
  {
    PRECICE_CHECK(math::greater(_r, 0.0),
                  "Support radius for radial-basis-function compact thin-plate-splines c2 has to be larger than zero. Please update the \"support-radius\" attribute.");
  }

  double getSupportRadius() const
  {
    return _r;
  }

  double evaluate(double radius) const
  {
    if (radius >= _r)
      return 0.0;
    double const p = radius / _r;
    using std::log;
    using std::pow;
    return 1.0 - 30.0 * pow(p, 2.0) - 10.0 * pow(p, 3.0) + 45.0 * pow(p, 4.0) - 6.0 * pow(p, 5.0) - 60.0 * log(pow(p, pow(p, 3.0)));
  }

private:
  logging::Logger _log{"mapping::CompactThinPlateSplinesC2"};

  double const _r;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^2,
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC0 : public CompactSupportBase {
public:
  explicit CompactPolynomialC0(double supportRadius)
      : _r(supportRadius)
  {
    PRECICE_CHECK(math::greater(_r, 0.0),
                  "Support radius for radial-basis-function compact polynomial c0 has to be larger than zero. Please update the \"support-radius\" attribute.");
  }

  double getSupportRadius() const
  {
    return _r;
  }

  double evaluate(double radius) const
  {
    if (radius >= _r)
      return 0.0;
    return std::pow(1.0 - radius / _r, 2.0);
  }

private:
  logging::Logger _log{"mapping::CompactPolynomialC0"};

  double const _r;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^8 * (32*rn^3 + 25*rn^2 + 8*rn + 1),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactPolynomialC6 : public CompactSupportBase {
public:
  explicit CompactPolynomialC6(double supportRadius)
      : _r(supportRadius)
  {
    PRECICE_CHECK(math::greater(_r, 0.0),
                  "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute.");
  }

  double getSupportRadius() const
  {
    return _r;
  }

  double evaluate(double radius) const
  {
    if (radius >= _r)
      return 0.0;
    double p = radius / _r;
    using std::pow;
    return pow(1.0 - p, 8.0) * (32.0 * pow(p, 3.0) + 25.0 * pow(p, 2.0) + 8.0 * p + 1.0);
  }

private:
  logging::Logger _log{"mapping::CompactPolynomialC6"};

  double const _r;
};

} // namespace mapping
} // namespace precice

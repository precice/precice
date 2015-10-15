#pragma once

#include "utils/Globals.hpp"

namespace precice {
namespace mapping {

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius^2 * log(radius).
 */
class ThinPlateSplines
{
public:

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    double result = 0.0;
    if (tarch::la::greater(radius, 0.0)){
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
 * Evaluates to: sqrt(shape + radius^2).
 */
class Multiquadrics
{
public:

  Multiquadrics ( double c )
    : _cPow2(std::pow(c, 2)) {}

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
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
class InverseMultiquadrics
{
public:

  InverseMultiquadrics ( double c )
    : _cPow2(std::pow(c, 2))
  {
    preciceCheck(tarch::la::greater(c, 0.0), "InverseMultiquadrics()",
                 "Shape parameter for radial-basis-function inverse multiquadric"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return 1.0 / std::sqrt(_cPow2 + std::pow(radius, 2));
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _cPow2;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius.
 */
class VolumeSplines
{
public:

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return radius;
  }
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: exp(-1 * shape * radius^2).
 */
class Gaussian
{
public:

  Gaussian ( double shape )
    : _shape(shape)
  {
    preciceCheck(tarch::la::greater(_shape, 0.0), "Gaussian()",
                 "Shape parameter for radial-basis-function gaussian"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return std::exp( - std::pow(_shape*radius,2.0) );
  }

private:

  static tarch::logging::Log _log;

  double _shape;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: 1 - 30*rn^2 - 10*rn^3 + 45*rn^4 - 6*rn^5 - 60*log(rn^3),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class CompactThinPlateSplinesC2
{
public:

  CompactThinPlateSplinesC2 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactThinPlateSplinesC2()",
                 "Support radius for radial-basis-function compact thin-plate-splines c2"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    double p = radius / _r;
    using std::pow;
    using std::log;
    return 1.0 - 30.0*pow(p,2.0) - 10.0*pow(p,3.0) + 45.0*pow(p,4.0)
      - 6.0*pow(p,5.0) - 60.0*log(pow(p,pow(p,3.0)));
  }

private:

  static tarch::logging::Log _log;

  double _r;
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
class CompactPolynomialC0
{
public:

  CompactPolynomialC0 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactPolynomialC0()",
                 "Support radius for radial-basis-function compact polynomial c0"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    return std::pow(1.0 - radius/_r, 2.0);
  }

private:

  static tarch::logging::Log _log;

  double _r;
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
class CompactPolynomialC6
{
public:

  CompactPolynomialC6 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactPolynomialC6()",
                 "Support radius for radial-basis-function compact polynomial c6"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    double p = radius / _r;
    using std::pow;
    return pow(1.0-p,8.0) * (32.0*pow(p,3.0) + 25.0*pow(p,2.0) + 8.0*p + 1.0);
  }

private:

  static tarch::logging::Log _log;

  double _r;
};

}} // namespace precice, mapping

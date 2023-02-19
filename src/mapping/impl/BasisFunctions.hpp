#pragma once

#include <array>
#include "logging/Logger.hpp"
#include "math/math.hpp"

#ifndef PRECICE_NO_GINKGO

#include "mapping/impl/DeviceBasisFunctions.cuh"

#endif

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
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const ThinPlateSplinesFunctor getFunctor() const
  {
    return this->_functor;
  }

private:
  std::array<double, 3>         _params;
  const std::string             _name = "ThinPlateSplines";
  const ThinPlateSplinesFunctor _functor{};

#endif
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
      : _cPow2(std::pow(c, 2))
  {
    _params.at(0) = _cPow2;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const MultiQuadraticsFunctor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  double                       _cPow2;
  const MultiQuadraticsFunctor _functor{};
  std::array<double, 3>        _params;
  const std::string            _name = "Multiquadratics";
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
      : _cPow2(std::pow(c, 2))
  {
    PRECICE_CHECK(math::greater(c, 0.0),
                  "Shape parameter for radial-basis-function inverse multiquadric has to be larger than zero. Please update the \"shape-parameter\" attribute.");
    _params.at(0) = _cPow2;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const InverseMultiquadricsFunctor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  logging::Logger _log{"mapping::InverseMultiQuadrics"};

  double const _cPow2;

  const InverseMultiquadricsFunctor _functor{};
  std::array<double, 3>             _params;
  const std::string                 _name = "InverseMultiquadratics";
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
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const VolumeSplinesFunctor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  const VolumeSplinesFunctor _functor{};
  std::array<double, 3>      _params;
  const std::string          _name = "VolumeSplines";
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
    PRECICE_CHECK(math::greater(_shape, 0.0),
                  "Shape parameter for radial-basis-function gaussian has to be larger than zero. Please update the \"shape-parameter\" attribute.");
    PRECICE_CHECK(math::greater(_supportRadius, 0.0),
                  "Support radius for radial-basis-function gaussian has to be larger than zero. Please update the \"support-radius\" attribute.");

    if (supportRadius < std::numeric_limits<double>::infinity()) {
      _deltaY = evaluate(supportRadius);
    }
    double threshold = std::sqrt(-std::log(cutoffThreshold)) / shape;
    _supportRadius   = std::min(supportRadius, threshold);
    _params.at(0)    = _shape;
    _params.at(1)    = _supportRadius;
    _params.at(2)    = _deltaY;
  }

  double getSupportRadius() const
  {
    return _supportRadius;
  }

  double evaluate(const double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  // TODO: Make precomputed
  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  };

  const std::string getName() const
  {
    return this->_name;
  }

  const GaussianFunctor getFunctor() const
  {
    return this->_functor;
  }

#endif
public:
  /// Below that value the function is supposed to be zero. Defines the support radius if not explicitly given
  static constexpr double cutoffThreshold = 1e-9;

private:
  logging::Logger _log{"mapping::Gaussian"};

  double const _shape;

  /// Either explicitly set (from cutoffThreshold) or computed supportRadius
  double _supportRadius;

  double _deltaY = 0;

  const GaussianFunctor _functor{};
  std::array<double, 3> _params;
  const std::string     _name = "Gaussian";
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
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact thin-plate-splines c2 has to be larger than zero. Please update the \"support-radius\" attribute.");
    _r_inv        = 1. / supportRadius;
    _params.at(0) = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const CompactThinPlateSplinesC2Functor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  logging::Logger _log{"mapping::CompactThinPlateSplinesC2"};

  double                                 _r_inv;
  const CompactThinPlateSplinesC2Functor _functor{};
  std::array<double, 3>                  _params;
  const std::string                      _name = "CompactThinPlateSplinesC2";
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
    logging::Logger _log{"mapping::CompactPolynomialC0"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c0 has to be larger than zero. Please update the \"support-radius\" attribute.");
    _r_inv        = 1. / supportRadius;
    _params.at(0) = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const CompactPolynomialC0Functor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  double                           _r_inv;
  const CompactPolynomialC0Functor _functor{};
  std::array<double, 3>            _params;
  const std::string                _name = "CompactPolynomialC0";
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
    logging::Logger _log{"mapping::CompactPolynomialC2"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c2 has to be larger than zero. Please update the \"support-radius\" attribute.");

    _r_inv        = 1. / supportRadius;
    _params.at(0) = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const CompactPolynomialC2Functor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  double                           _r_inv;
  const CompactPolynomialC2Functor _functor{};
  std::array<double, 3>            _params;
  const std::string                _name = "CompactPolynomialC2";
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
    logging::Logger _log{"mapping::CompactPolynomialC4"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c4 has to be larger than zero. Please update the \"support-radius\" attribute.");

    _r_inv        = 1. / supportRadius;
    _params.at(0) = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  std::array<double, 3> getFunctionParameters()
  {
    return _params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const CompactPolynomialC4Functor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  double                           _r_inv;
  const CompactPolynomialC4Functor _functor{};
  std::array<double, 3>            _params;
  const std::string                _name = "CompactPolynomialC4";
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
    logging::Logger _log{"mapping::CompactPolynomialC6"};
    PRECICE_CHECK(math::greater(supportRadius, 0.0),
                  "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute.");

    _r_inv              = 1. / supportRadius;
    this->_params.at(0) = _r_inv;
  }

  double getSupportRadius() const
  {
    return 1. / _r_inv;
  }

  double evaluate(double radius) const
  {
    return this->_functor(radius, this->_params);
  }

#ifndef PRECICE_NO_GINKGO

  const std::array<double, 3> getFunctionParameters()
  {
    return this->_params;
  }

  const std::string getName() const
  {
    return this->_name;
  }

  const CompactPolynomialC6Functor getFunctor() const
  {
    return this->_functor;
  }

#endif

private:
  double                           _r_inv;
  const CompactPolynomialC6Functor _functor{};
  std::array<double, 3>            _params;
  const std::string                _name = "CompactPolynomialC6";
};
} // namespace mapping
} // namespace precice

#include "mapping/impl/BasisFunctions.hpp"

#include "logging/Logger.hpp"
#include "math/math.hpp"

namespace precice {
namespace mapping {

InverseMultiquadrics ::InverseMultiquadrics(double c)
    : _cPow2(math::pow_int<2>(c))
{
  logging::Logger _log{"mapping::InverseMultiQuadrics"};
  PRECICE_CHECK(math::greater(c, 0.0),
                "Shape parameter for radial-basis-function inverse multiquadric has to be larger than zero. Please update the \"shape-parameter\" attribute.");
  _params.parameter1 = _cPow2;
}

Gaussian::Gaussian(const double shape, const double supportRadius)
    : _shape(shape),
      _supportRadius(supportRadius)
{
  logging::Logger _log{"mapping::Gaussian"};
  PRECICE_CHECK(math::greater(_shape, 0.0),
                "Shape parameter for radial-basis-function gaussian has to be larger than zero. Please update the \"shape-parameter\" attribute.");
  PRECICE_CHECK(math::greater(_supportRadius, 0.0),
                "Support radius for radial-basis-function gaussian has to be larger than zero. Please update the \"support-radius\" attribute.");

  double threshold   = std::sqrt(-std::log(cutoffThreshold)) / shape;
  _supportRadius     = std::min(supportRadius, threshold);
  _params.parameter1 = _shape;
  _params.parameter2 = _supportRadius;
  _params.parameter3 = _deltaY;
  if (supportRadius < std::numeric_limits<double>::infinity()) {
    _deltaY            = evaluate(supportRadius);
    _params.parameter3 = _deltaY;
  }
}

CompactThinPlateSplinesC2::CompactThinPlateSplinesC2(double supportRadius)
{
  logging::Logger _log{"mapping::CompactThinPlateSplinesC2"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact thin-plate-splines c2 has to be larger than zero. Please update the \"support-radius\" attribute.");
  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

CompactPolynomialC0 ::CompactPolynomialC0(double supportRadius)
{
  logging::Logger _log{"mapping::CompactPolynomialC0"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact polynomial c0 has to be larger than zero. Please update the \"support-radius\" attribute.");
  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

CompactPolynomialC2 ::CompactPolynomialC2(double supportRadius)
{
  logging::Logger _log{"mapping::CompactPolynomialC2"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact polynomial c2 has to be larger than zero. Please update the \"support-radius\" attribute.");

  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

CompactPolynomialC4 ::CompactPolynomialC4(double supportRadius)
{
  logging::Logger _log{"mapping::CompactPolynomialC4"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact polynomial c4 has to be larger than zero. Please update the \"support-radius\" attribute.");

  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

CompactPolynomialC6::CompactPolynomialC6(double supportRadius)
{
  logging::Logger _log{"mapping::CompactPolynomialC6"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute.");
  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

CompactPolynomialC8 ::CompactPolynomialC8(double supportRadius)
{
  logging::Logger _log{"mapping::CompactPolynomialC8"};
  PRECICE_CHECK(math::greater(supportRadius, 0.0),
                "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute.");
  _r_inv             = 1. / supportRadius;
  _params.parameter1 = _r_inv;
}

} // namespace mapping
} // namespace precice

#ifndef PRECICE_NO_GINKGO

#include <cmath>
#include <stdio.h>

#include "mapping/impl/DeviceBasisFunctions.cuh"

#define NUMERICAL_ZERO_DIFFERENCE 1.0e-14

namespace precice {
namespace mapping {

SHARED_HOST_DEVICE_FUNCTION double ThinPlateSplinesFunctor::operator()(const double radius, const std::array<double, 3> params) const
{
  // We don't need to read any values from params since there is no need here
  return std::log(std::max(radius, NUMERICAL_ZERO_DIFFERENCE)) * std::pow(radius, 2);
}

SHARED_HOST_DEVICE_FUNCTION double MultiQuadraticsFunctor::operator()(const double radius, const std::array<double, 3> params) const
{
  double cPow2 = params.at(0);
  return std::sqrt(cPow2 + std::pow(radius, 2));
}

SHARED_HOST_DEVICE_FUNCTION double InverseMultiquadricsFunctor::operator()(const double radius, const std::array<double, 3> params) const
{
  double cPow2 = params.at(0);
  return 1.0 / std::sqrt(cPow2 + std::pow(radius, 2));
}

SHARED_HOST_DEVICE_FUNCTION double VolumeSplinesFunctor::operator()(const double radius, const std::array<double, 3> params) const
{
  return std::abs(radius);
}

SHARED_HOST_DEVICE_FUNCTION double GaussianFunctor::operator()(const double radius, std::array<double, 3> params) const
{
  double shape         = params.at(0);
  double supportRadius = params.at(1);
  double deltaY        = params.at(2);

  if (radius > supportRadius) {
    return 0.0;
  } else {
    return std::exp(-std::pow(shape * radius, 2)) - deltaY;
  }
}

SHARED_HOST_DEVICE_FUNCTION double CompactThinPlateSplinesC2Functor::operator()(const double radius, const std::array<double, 3> params) const
{
  double       r_inv = params.at(0);
  double const p     = radius * r_inv;
  if (p >= 1) {
    return 0.0;
  } else {
    return 1.0 - 30.0 * std::pow(p, 2) - 10.0 * std::pow(p, 3) + 45.0 * std::pow(p, 4) - 6.0 * std::pow(p, 5) - std::pow(p, 3) * 60.0 * std::log(std::max(p, NUMERICAL_ZERO_DIFFERENCE));
  }
}

SHARED_HOST_DEVICE_FUNCTION double CompactPolynomialC0Functor::operator()(const double radius, const std::array<double, 3> params) const
{
  double       r_inv = params.at(0);
  double const p     = radius * r_inv;
  if (p >= 1) {
    return 0.0;
  } else {
    return std::pow(1.0 - p, 2);
  }
}

SHARED_HOST_DEVICE_FUNCTION double CompactPolynomialC2Functor::operator()(const double radius, const std::array<double, 3> params) const
{
  double       r_inv = params.at(0);
  double const p     = radius * r_inv;
  if (p >= 1) {
    return 0.0;
  } else {
    return std::pow(1.0 - p, 4) * (4 * p + 1);
  }
}

SHARED_HOST_DEVICE_FUNCTION double CompactPolynomialC4Functor::operator()(const double radius, const std::array<double, 3> params) const
{
  double       r_inv = params.at(0);
  double const p     = radius * r_inv;
  if (p >= 1) {
    return 0.0;
  } else {
    return std::pow(1.0 - p, 6) * (35 * std::pow(p, 2) + 18 * p + 3, 2);
  }
}

SHARED_HOST_DEVICE_FUNCTION double CompactPolynomialC6Functor::operator()(const double radius, const std::array<double, 3> params) const
{
  double       r_inv = params.at(0);
  double const p     = radius * r_inv;
  if (p >= 1) {
    return 0.0;
  } else {
    return std::pow(1.0 - p, 8) * (32.0 * std::pow(p, 3) + 25.0 * std::pow(p, 2) + 8.0 * p + 1.0);
  }
}

} // namespace mapping
} // namespace precice

#endif

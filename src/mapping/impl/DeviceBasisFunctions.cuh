#pragma once
#ifndef PRECICE_NO_GINKGO

#ifdef __NVCC__

#include <cuda_runtime.h>
#define SHARED_HOST_DEVICE_FUNCTION __host__ __device__

#else

#define SHARED_HOST_DEVICE_FUNCTION

#endif

#include <array>
#include <string>

namespace precice {
namespace mapping {

class ThinPlateSplinesFunctor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class MultiQuadraticsFunctor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class InverseMultiquadricsFunctor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class VolumeSplinesFunctor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class GaussianFunctor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class CompactThinPlateSplinesC2Functor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class CompactPolynomialC0Functor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class CompactPolynomialC2Functor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class CompactPolynomialC4Functor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

class CompactPolynomialC6Functor {
public:
  SHARED_HOST_DEVICE_FUNCTION double operator()(const double radius, const std::array<double, 3> params) const;
};

} // namespace mapping
} // namespace precice

#endif

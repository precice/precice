#pragma once

/**
 * @file BasisFunctions.hpp
 * @brief Radial Basis Function (RBF) kernel definitions for data mapping in preCICE.
 *
 * ## What are Radial Basis Functions?
 * An RBF kernel phi(r) is a function that depends only on the Euclidean distance r = ||x - x_i||
 * between two points x and x_i. Given input data values f(x_i) at a set of N input vertices,
 * the RBF interpolant is:
 *
 *   s(x) = sum_{i=1}^{N} lambda_i * phi(||x - x_i||) + p(x)
 *
 * where:
 *   - lambda_i are the interpolation coefficients (solved from a linear system)
 *   - phi is the chosen radial basis function (defined in this file)
 *   - p(x) is a low-order polynomial (constant + linear terms for conditionally positive definite RBFs)
 *
 * ## Choosing a Basis Function
 * - **Compact-support RBFs** (Gaussian, Wendland C0/C2/C4/C6/C8, CompactTPS-C2):
 *   Only interact with nearby vertices ⟹ sparse RBF matrix ⟹ faster for large meshes.
 *   Require a `support-radius` attribute in the preCICE config.
 *
 * - **Global-support RBFs** (ThinPlateSplines, VolumeSplines, Multiquadrics, InverseMultiquadrics):
 *   Interact with ALL vertices ⟹ dense RBF matrix ⟹ accurate but slow for large meshes.
 *   No support radius needed.
 *
 * - **Strictly positive definite RBFs** (InverseMultiquadrics, Gaussian, Wendland families):
 *   The RBF matrix C is symmetric positive definite ⟹ efficient Cholesky decomposition can be used.
 *   These functions do NOT support the integrated polynomial option (polynomial="on").
 *
 * - **Conditionally positive definite RBFs** (ThinPlateSplines, VolumeSplines, Multiquadrics):
 *   Require an additional polynomial for well-posedness.
 *   The RBF matrix is decomposed via QR factorization (more general, slightly slower).
 *
 * ## GPU Compatibility
 * All basis function evaluation operators are annotated with `PRECICE_HOST_DEVICE`,
 * making them callable both on the CPU and on GPU device kernels (CUDA/HIP).
 * The platform-specific intrinsics (fma, log) are abstracted via PRECICE_FMA and PRECICE_LOG macros.
 */

#include <stdexcept>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

#ifdef __HIPCC__
#include <hip/hip_runtime.h>
#endif

#if defined(__CUDACC__) || defined(__HIPCC__)

#include <ginkgo/extensions/kokkos.hpp>
#include <ginkgo/ginkgo.hpp>

#define PRECICE_HOST_DEVICE __host__ __device__
#define PRECICE_MEMORY_SPACE __device__
#define PRECICE_FMA Kokkos::fma
#define PRECICE_LOG Kokkos::log

#else

#define PRECICE_HOST_DEVICE
#define PRECICE_MEMORY_SPACE
#define PRECICE_FMA std::fma
#define PRECICE_LOG std::log

#endif

#define NUMERICAL_ZERO_DIFFERENCE_DEVICE 1.0e-14

#include "math/math.hpp"

namespace precice::mapping {

/**
 * @brief Wrapper struct that is used to transfer RBF-specific parameters to the GPU.
 *
 * Its parameters 1, 2 and 3 are filled individually, depending on the requirements
 * of each RBF (e.g., support radius, shape parameter, etc.).
 * This parameter struct is handed to the GPU device kernel
 * when assembling the interpolation matrices
 * and its parameters used in every RBF evaluation.
 * Since there are at most 3 parameters for a RBF in the
 * current implementation, 3 parameters are required in this
 * struct, but ignored if the specific RBF does not require them all.
 *
 */
struct RadialBasisParameters {
  /// parameter1: shape parameter squared (Multiquadrics, InverseMultiquadrics),
  ///             shape parameter (Gaussian), or 1/supportRadius (compact Wendland/TPS families).
  double parameter1{};
  /// parameter2: support radius (Gaussian only — stores effective cutoff radius after threshold).
  double parameter2{};
  /// parameter3: deltaY shift (Gaussian only — vertical offset to ensure phi(supportRadius) → 0 smoothly).
  double parameter3{};
};

/**
 * @brief Trait base class: marks an RBF as having compact (local) support.
 *
 * Compact-support RBFs evaluate to exactly zero beyond a finite support radius.
 * This produces a sparse interpolation matrix, which is significantly faster for large meshes.
 * Derived classes must also provide a `getSupportRadius()` method.
 */
struct CompactSupportBase {
  static constexpr bool hasCompactSupport()
  {
    return true;
  }
};

/**
 * @brief Trait base class: marks an RBF as having global (unbounded) support.
 *
 * Global-support RBFs are non-zero everywhere in the domain.
 * This leads to a fully dense kernel matrix — expensive for large meshes,
 * but can be more accurate for small-to-medium problems.
 * getSupportRadius() returns the maximum double value to signal "no cutoff".
 */
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
 * @brief Trait base class: encodes whether the RBF is strictly positive definite.
 *
 * @tparam isDefinite true if the kernel produces a symmetric positive definite (SPD) matrix.
 *
 * This distinction determines which matrix decomposition is used in RadialBasisFctSolver:
 * - `isDefinite = true`  → Cholesky decomposition (LLT): fast, numerically stable for SPD matrices.
 * - `isDefinite = false` → QR decomposition (ColPivHouseholderQR): works for non-SPD, requires polynomial.
 *
 * Strictly positive definite kernels (InverseMultiquadrics, Gaussian, Wendland families)
 * do NOT support the integrated polynomial mode (polynomial="on"), as the augmented system
 * loses its positive definiteness structure.
 */
template <bool isDefinite>
struct DefiniteFunction {
  static constexpr bool isStrictlyPositiveDefinite()
  {
    return isDefinite;
  }
};

/**
 * @brief Thin-Plate Spline (TPS) — globally supported, conditionally positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula:  phi(r) = r^2 * log(r)
 *
 * Properties:
 * - Global support: interacts with ALL vertices → dense matrix.
 * - Conditionally positive definite (order 2) → requires a polynomial (use polynomial="separate" or "on").
 * - No shape/support-radius parameter needed.
 * - Smooth (C^inf away from r=0), well-suited for smooth fields.
 *
 * Note: log(0) is handled by clamping r to NUMERICAL_ZERO_DIFFERENCE_DEVICE before evaluation.
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
    return PRECICE_LOG(std::max(radius, NUMERICAL_ZERO_DIFFERENCE_DEVICE)) * math::pow_int<2>(radius);
  }

  RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  RadialBasisParameters _params;
};

/**
 * @brief Multiquadric — globally supported, conditionally positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula:  phi(r) = sqrt(c^2 + r^2),   c = shape parameter > 0
 *
 * Properties:
 * - Global support → dense matrix.
 * - Conditionally positive definite → requires polynomial (use polynomial="separate" or "on").
 * - Config attribute: `shape-parameter` (stored as c^2 in parameter1).
 * - Larger c → smoother interpolant; smaller c → more local behaviour.
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
 * @brief Inverse Multiquadric — globally supported, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula:  phi(r) = 1 / sqrt(c^2 + r^2),   c = shape parameter > 0
 *
 * Properties:
 * - Global support → dense matrix.
 * - Strictly positive definite → Cholesky decomposition can be used (efficient).
 * - Config attribute: `shape-parameter` (stored as c^2 in parameter1).
 * - Does NOT support polynomial="on" (integrated polynomial) due to SPD requirement.
 * - Decays to 0 as r → ∞, unlike Multiquadrics which grows.
 */
class InverseMultiquadrics : public NoCompactSupportBase,
                             public DefiniteFunction<true> {
public:
  explicit InverseMultiquadrics(double c)
      : _cPow2(math::pow_int<2>(c))
  {
    if (!math::greater(c, 0.0)) {
      throw std::invalid_argument{"Shape parameter for radial-basis-function inverse multiquadric has to be larger than zero. Please update the \"shape-parameter\" attribute."};
    }
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
 * @brief Volume Spline — globally supported, conditionally positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula:  phi(r) = |r|
 *
 * Properties:
 * - Global support → dense matrix.
 * - Conditionally positive definite (order 1) → requires polynomial.
 * - No shape/support-radius parameter needed.
 * - Simpler than TPS but less smooth (only C^1); useful for piecewise-linear-like behaviour.
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
 * @brief Gaussian — can be compact or global, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula:  phi(r) = exp(-(shape * r)^2) - deltaY
 *   where deltaY = phi(supportRadius) ensures the function smoothly reaches 0 at the boundary.
 *
 * Config attributes:
 *   - `shape-parameter` (shape > 0): controls the width of the Gaussian bell.
 *     Large shape → narrow bell (local behaviour); small shape → wide bell (global-like).
 *   - `support-radius` (optional): if given, the Gaussian is truncated at this radius.
 *     If not given, the cutoff threshold (1e-9) is used to automatically determine the radius.
 *
 * Properties:
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - parameter1 = shape, parameter2 = effective support radius, parameter3 = deltaY shift.
 * - Very smooth (C^inf) but can be sensitive to the shape parameter choice.
 */
class Gaussian : public CompactSupportBase,
                 public DefiniteFunction<true> {
public:
  Gaussian(const double shape, const double supportRadius = std::numeric_limits<double>::infinity())
      : _shape(shape),
        _supportRadius(supportRadius)
  {
    if (!math::greater(_shape, 0.0)) {
      throw std::invalid_argument{
          "Shape parameter for radial-basis-function gaussian has to be larger than zero. Please update the \"shape-parameter\" attribute."};
    }
    if (!math::greater(_supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function gaussian has to be larger than zero. Please update the \"support-radius\" attribute."};
    }

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
 * @brief Compact Thin-Plate Spline C2 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only evaluated for rn < 1):
 *   phi(r) = 1 - 30*rn^2 - 10*rn^3 + 45*rn^4 - 6*rn^5 - 60*rn^3 * log(rn)
 *
 * The log(rn) term is carefully handled to avoid log(0) when r=0.
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius = _r_inv in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix for large meshes.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^2 continuity: second derivatives exist everywhere including at r=0.
 */
class CompactThinPlateSplinesC2 : public CompactSupportBase,
                                  public DefiniteFunction<true> {
public:
  explicit CompactThinPlateSplinesC2(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact thin-plate-splines c2 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }
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
    return 1.0 - 30.0 * math::pow_int<2>(p) - 10.0 * math::pow_int<3>(p) + 45.0 * math::pow_int<4>(p) - 6.0 * math::pow_int<5>(p) - math::pow_int<3>(p) * 60.0 * PRECICE_LOG(std::max(p, NUMERICAL_ZERO_DIFFERENCE_DEVICE));
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
 * @brief Wendland C0 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only for rn < 1):
 *   phi(r) = (1 - rn)^2
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^0 continuity: continuous but NOT differentiable at r=0.
 *   Suitable for non-smooth fields or for a simpler/cheaper kernel.
 */
class CompactPolynomialC0 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC0(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact polynomial c0 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }
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
 * @brief Wendland C2 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only for rn < 1):
 *   phi(r) = (1 - rn)^4 * (4*rn + 1)
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^2 continuity: smooth up to second derivatives everywhere.
 *   Good default choice for smooth physical fields (temperature, pressure, displacement).
 */
class CompactPolynomialC2 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC2(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact polynomial c2 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }

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
    return math::pow_int<4>(1.0 - p) * PRECICE_FMA(4, p, 1);
  }

  RadialBasisParameters getFunctionParameters() const
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland C4 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only for rn < 1):
 *   phi(r) = (1 - rn)^6 * (35*rn^2 + 18*rn + 3)
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^4 continuity: four times continuously differentiable.
 *   Smoother than C2; use when interpolating highly smooth fields.
 */
class CompactPolynomialC4 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC4(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact polynomial c4 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }

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
    return math::pow_int<6>(1.0 - p) * (35 * math::pow_int<2>(p) + PRECICE_FMA(18, p, 3));
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
 * @brief Wendland C6 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only for rn < 1):
 *   phi(r) = (1 - rn)^8 * (32*rn^3 + 25*rn^2 + 8*rn + 1)
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^6 continuity: six times continuously differentiable.
 *   Very smooth; use for highly smooth fields when maximal regularity is needed.
 */
class CompactPolynomialC6 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC6(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }
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
    return math::pow_int<8>(1.0 - p) * (32.0 * math::pow_int<3>(p) + 25.0 * math::pow_int<2>(p) + PRECICE_FMA(8.0, p, 1.0));
  };

  const RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

/**
 * @brief Wendland C8 — compact support, strictly positive definite.
 *
 * Use as a template parameter for RadialBasisFctMapping.
 *
 * Formula (rn = r / supportRadius, only for rn < 1):
 *   phi(r) = (1 - rn)^10 * (1287*rn^4 + 1350*rn^3 + 630*rn^2 + 150*rn + 15)
 *
 * Config attribute: `support-radius` (stored as 1/supportRadius in parameter1).
 *
 * Properties:
 * - Compact support → sparse matrix.
 * - Strictly positive definite → Cholesky decomposition. Does NOT support polynomial="on".
 * - C^8 continuity: eight times continuously differentiable.
 *   The smoothest of the Wendland family implemented here.
 *   Use for very smooth fields or when high-order accuracy is critical.
 */
class CompactPolynomialC8 : public CompactSupportBase,
                            public DefiniteFunction<true> {
public:
  explicit CompactPolynomialC8(double supportRadius)
  {
    if (!math::greater(supportRadius, 0.0)) {
      throw std::invalid_argument{
          "Support radius for radial-basis-function compact polynomial c6 has to be larger than zero. Please update the \"support-radius\" attribute."};
    }
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
    return math::pow_int<10>(1.0 - p) * (1287.0 * math::pow_int<4>(p) + 1350.0 * math::pow_int<3>(p) + 630.0 * math::pow_int<2>(p) + 150.0 * p + 15);
  };

  const RadialBasisParameters getFunctionParameters()
  {
    return _params;
  }

private:
  double                _r_inv;
  RadialBasisParameters _params;
};

#undef PRECICE_FMA
#undef PRECICE_LOG
#undef PRECICE_MEMORY_SPACE
#undef PRECICE_HOST_DEVICE
#undef NUMERICAL_ZERO_DIFFERENCE_DEVICE

} // namespace precice::mapping

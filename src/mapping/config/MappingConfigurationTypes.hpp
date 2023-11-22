#pragma once

namespace precice::mapping {

/// How to handle the polynomial?
/**
 * ON: Include it in the system matrix
 * OFF: Omit it altogether
 * SEPARATE: Compute it separately using least-squares QR.
 */
enum class Polynomial {
  ON,
  OFF,
  SEPARATE
};

enum class BasisFunction {
  WendlandC0,
  WendlandC2,
  WendlandC4,
  WendlandC6,
  WendlandC8,
  ThinPlateSplines,
  Multiquadrics,
  InverseMultiquadrics,
  VolumeSplines,
  Gaussian,
  CompactThinPlateSplinesC2,
};
} // namespace precice::mapping

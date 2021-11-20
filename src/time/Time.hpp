#pragma once

namespace precice {
namespace time {

class Time {
public:
  /// To be used, when the extrapolation order is not defined.
  static const int UNDEFINED_EXTRAPOLATION_ORDER;

  /// To be used, when the interpolation order is not defined.
  static const int UNDEFINED_INTERPOLATION_ORDER;
};

} // namespace time
} // namespace precice
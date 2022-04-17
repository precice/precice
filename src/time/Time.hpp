#pragma once

namespace precice {
namespace time {

class Time {
public:
  /// To be used, when the interpolation order is not defined.
  static const int DEFAULT_INTERPOLATION_ORDER;

  /// The minimum required interpolation order.
  static const int MIN_INTERPOLATION_ORDER;

  /// The maximum allowed interpolation order.
  static const int MAX_INTERPOLATION_ORDER;
};

} // namespace time
} // namespace precice

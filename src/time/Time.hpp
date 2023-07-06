#pragma once

namespace precice {
namespace time {

class Time {
public:
  /// To be used, when the interpolation order is not defined.
  static const int DEFAULT_WAVEFORM_DEGREE;

  /// The minimum required interpolation order.
  static const int MIN_WAVEFORM_DEGREE;

  /// The maximum allowed interpolation order.
  static const int MAX_WAVEFORM_DEGREE;
};

} // namespace time
} // namespace precice

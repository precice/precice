#pragma once

namespace precice {
namespace time {

class Time {
public:
  /// To be used, when the interpolation degree is not defined.
  static const int DEFAULT_WAVEFORM_DEGREE;

  /// The minimum required interpolation degree.
  static const int MIN_WAVEFORM_DEGREE;
};

} // namespace time
} // namespace precice

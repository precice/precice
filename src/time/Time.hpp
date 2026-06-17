#pragma once

namespace precice::time {

class Time {
public:
  /// To be used, when the interpolation degree is not defined.
  inline static const int DEFAULT_WAVEFORM_DEGREE = 1;

  /// The minimum required interpolation degree.
  inline static const int MIN_WAVEFORM_DEGREE = 0;
};

} // namespace precice::time

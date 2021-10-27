#pragma once

#include "time/Waveform.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the Waveform class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
class WaveformFixture {
public:
  int numberOfStoredSamples(time::Waveform &waveform);

  int sizeOfSampleStorage(time::Waveform &waveform);

  int dataCount(time::Waveform &waveform);
};

} // namespace testing
} // namespace precice

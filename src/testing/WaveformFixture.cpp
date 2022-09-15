#include "testing/WaveformFixture.hpp"

namespace precice::testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._timeStepsStorage.size();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform._timeStepsStorage[0.0].size();
}

} // namespace precice::testing

#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._timeStepsStorage.nTimes();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform._timeStepsStorage.nDofs();
}

} // namespace testing
} // namespace precice

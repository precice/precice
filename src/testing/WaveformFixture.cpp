#include "testing/WaveformFixture.hpp"
#include "mesh/Data.hpp"

namespace precice::testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._data->timeStepsStorage().nTimes();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform._data->timeStepsStorage().nDofs();
}

} // namespace precice::testing

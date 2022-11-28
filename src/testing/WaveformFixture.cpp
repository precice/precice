#include "testing/WaveformFixture.hpp"

namespace precice::testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._storage.nTimes();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform._storage.nDofs();
}

} // namespace precice::testing

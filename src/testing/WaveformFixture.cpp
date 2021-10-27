#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._numberOfStoredSamples;
}

int WaveformFixture::sizeOfSampleStorage(time::Waveform &waveform)
{
  return waveform.sizeOfSampleStorage();
}

int WaveformFixture::dataCount(time::Waveform &waveform)
{
  return waveform.dataCount();
}
} // namespace testing
} // namespace precice

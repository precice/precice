#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfSamples(time::Waveform &waveform)
{
  return waveform.numberOfSamples();
}

int WaveformFixture::numberOfData(time::Waveform &waveform)
{
  return waveform.numberOfData();
}
} // namespace testing
} // namespace precice
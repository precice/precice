#include "testing/SerialCouplingSchemeFixture.hpp"

namespace precice {
namespace testing {

bool SerialCouplingSchemeFixture::isImplicitCouplingScheme(cplscheme::SerialCouplingScheme &cplscheme)
{
  return cplscheme.isImplicitCouplingScheme();
}

cplscheme::CouplingData *SerialCouplingSchemeFixture::getReceiveData(cplscheme::SerialCouplingScheme &cplscheme, DataID dataID)
{
  return cplscheme.getReceiveData(dataID);
}

cplscheme::CouplingData *SerialCouplingSchemeFixture::getSendData(cplscheme::SerialCouplingScheme &cplscheme, DataID dataID)
{
  return cplscheme.getSendData(dataID);
}

void SerialCouplingSchemeFixture::setTimeWindows(cplscheme::SerialCouplingScheme &cplscheme, int timeWindows)
{
  cplscheme.setTimeWindows(timeWindows);
}

void SerialCouplingSchemeFixture::storeIteration(cplscheme::SerialCouplingScheme &cplscheme)
{
  cplscheme.storeIteration();
}

void SerialCouplingSchemeFixture::initializeStorage(cplscheme::SerialCouplingScheme &cplscheme)
{
  cplscheme.initializeStorage();
}

void SerialCouplingSchemeFixture::storeDataInWaveforms(cplscheme::SerialCouplingScheme &cplscheme)
{
  cplscheme.storeDataInWaveforms();
}

void SerialCouplingSchemeFixture::moveToNextWindow(cplscheme::SerialCouplingScheme &cplscheme)
{
  cplscheme.moveToNextWindow();
}
} // namespace testing
} // namespace precice
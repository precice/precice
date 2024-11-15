#include "testing/SerialCouplingSchemeFixture.hpp"
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"

namespace precice::testing {

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

void SerialCouplingSchemeFixture::initializeAcceleration(cplscheme::SerialCouplingScheme &cplscheme)
{
  if (cplscheme._acceleration) {
    cplscheme._acceleration->initialize(cplscheme.getAccelerationData());
  }
}

void SerialCouplingSchemeFixture::moveToNextWindow(cplscheme::SerialCouplingScheme &cplscheme)
{
  for (const auto &pair : cplscheme._allData) {
    pair.second->timeStepsStorage().move();
  }
}
} // namespace precice::testing

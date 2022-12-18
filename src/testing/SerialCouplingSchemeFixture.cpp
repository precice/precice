#include "testing/SerialCouplingSchemeFixture.hpp"

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

void SerialCouplingSchemeFixture::initializeStorages(cplscheme::SerialCouplingScheme &cplscheme)
{
  cplscheme.initializeStorages();
}

void SerialCouplingSchemeFixture::moveToNextWindow(cplscheme::SerialCouplingScheme &cplscheme)
{
  for (const cplscheme::SerialCouplingScheme::DataMap::value_type &pair : cplscheme.getAllData()) {
    pair.second->moveTimeStepsStorage();
  }
}
} // namespace precice::testing

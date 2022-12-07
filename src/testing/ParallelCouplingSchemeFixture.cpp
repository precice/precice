#include "testing/ParallelCouplingSchemeFixture.hpp"

namespace precice::testing {

bool ParallelCouplingSchemeFixture::isImplicitCouplingScheme(cplscheme::ParallelCouplingScheme &cplscheme)
{
  return cplscheme.isImplicitCouplingScheme();
}

cplscheme::CouplingData *ParallelCouplingSchemeFixture::getReceiveData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getReceiveData(dataID);
}

cplscheme::CouplingData *ParallelCouplingSchemeFixture::getSendData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getSendData(dataID);
}

void ParallelCouplingSchemeFixture::setTimeWindows(cplscheme::ParallelCouplingScheme &cplscheme, int timeWindows)
{
  cplscheme.setTimeWindows(timeWindows);
}

void ParallelCouplingSchemeFixture::storeIteration(cplscheme::ParallelCouplingScheme &cplscheme)
{
  cplscheme.storeIteration();
}

void ParallelCouplingSchemeFixture::initializeStorages(cplscheme::ParallelCouplingScheme &cplscheme)
{
  cplscheme.initializeStorages();
}

void ParallelCouplingSchemeFixture::storeExtrapolationData(cplscheme::ParallelCouplingScheme &cplscheme)
{
  cplscheme.storeExtrapolationData();
}

void ParallelCouplingSchemeFixture::moveToNextWindow(cplscheme::ParallelCouplingScheme &cplscheme)
{
  cplscheme.moveToNextWindow();
}
} // namespace precice::testing

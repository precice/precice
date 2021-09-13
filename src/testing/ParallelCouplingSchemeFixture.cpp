#include "testing/ParallelCouplingSchemeFixture.hpp"

namespace precice {
namespace testing {

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
} // namespace testing
} // namespace precice
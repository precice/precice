#include "testing/SerialCouplingSchemeFixture.hpp"

namespace precice {
namespace testing {

bool SerialCouplingSchemeFixture::isImplicitCouplingScheme(cplscheme::SerialCouplingScheme &cplscheme)
{
  return cplscheme.isImplicitCouplingScheme();
}

cplscheme::CouplingData *SerialCouplingSchemeFixture::getReceiveData(cplscheme::SerialCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getReceiveData(dataID);
}

cplscheme::CouplingData *SerialCouplingSchemeFixture::getSendData(cplscheme::SerialCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getSendData(dataID);
}
} // namespace testing
} // namespace precice
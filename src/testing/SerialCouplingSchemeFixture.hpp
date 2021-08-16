#pragma once

#include "cplscheme/SerialCouplingScheme.hpp"

namespace precice {
    /*
namespace cplscheme {
// Forward declaration
class SerialCouplingScheme;
} // namespace cplscheme
*/
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the SerialCouplingScheme class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
class SerialCouplingSchemeFixture {
    public:
  bool isImplicitCouplingScheme(cplscheme::SerialCouplingScheme &cplscheme);

  cplscheme::CouplingData *getReceiveData(cplscheme::SerialCouplingScheme &cplscheme, int dataID);

  cplscheme::CouplingData *getSendData(cplscheme::SerialCouplingScheme &cplscheme, int dataID);
};

} // namespace testing
} // namespace precice
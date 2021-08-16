#pragma once

#include "cplscheme/ParallelCouplingScheme.hpp"

namespace precice {
/*
namespace cplscheme {
// Forward declaration
class ParallelCouplingScheme;
} // namespace cplscheme
*/
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the ParallelCouplingScheme class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
class ParallelCouplingSchemeFixture {
    public:
  bool isImplicitCouplingScheme(cplscheme::ParallelCouplingScheme &cplscheme);

  cplscheme::CouplingData *getReceiveData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID);

  cplscheme::CouplingData *getSendData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID);
};
} // namespace testing
} // namespace precice
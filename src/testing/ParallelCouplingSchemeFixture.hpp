#pragma once

#include "cplscheme/ParallelCouplingScheme.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the ParallelCouplingScheme class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
struct ParallelCouplingSchemeFixture {
  static bool isImplicitCouplingScheme(cplscheme::ParallelCouplingScheme &cplscheme);

  static cplscheme::CouplingData *getReceiveData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID);

  static cplscheme::CouplingData *getSendData(cplscheme::ParallelCouplingScheme &cplscheme, int dataID);
};
} // namespace testing
} // namespace precice

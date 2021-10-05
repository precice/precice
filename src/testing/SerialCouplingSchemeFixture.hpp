#pragma once

#include "cplscheme/SerialCouplingScheme.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the SerialCouplingScheme class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
struct SerialCouplingSchemeFixture {
  bool isImplicitCouplingScheme(cplscheme::SerialCouplingScheme &cplscheme);

  cplscheme::CouplingData *getReceiveData(cplscheme::SerialCouplingScheme &cplscheme, DataID dataID);

  cplscheme::CouplingData *getSendData(cplscheme::SerialCouplingScheme &cplscheme, DataID dataID);

  void setTimeWindows(cplscheme::SerialCouplingScheme &cplscheme, int timeWindows);

  void storeIteration(cplscheme::SerialCouplingScheme &cplscheme);

  void setupDataMatrices(cplscheme::SerialCouplingScheme &cplscheme);

  void storeDataInWaveforms(cplscheme::SerialCouplingScheme &cplscheme);

  void moveToNextWindow(cplscheme::SerialCouplingScheme &cplscheme);
};

} // namespace testing
} // namespace precice

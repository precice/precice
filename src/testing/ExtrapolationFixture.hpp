#pragma once

#include "cplscheme/impl/Extrapolation.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the Extrapolation class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
class ExtrapolationFixture {
public:
  int numberOfStoredSamples(cplscheme::Extrapolation &extrapolation);

  int sizeOfSampleStorage(cplscheme::Extrapolation &extrapolation);

  int valuesSize(cplscheme::Extrapolation &extrapolation);

  double getValue(cplscheme::Extrapolation &extrapolation, int valueID, int sampleID);
};

} // namespace testing
} // namespace precice

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
  int numberOfStoredSamples(cplscheme::impl::Extrapolation &extrapolation);

  int sizeOfSampleStorage(cplscheme::impl::Extrapolation &extrapolation);

  int valuesSize(cplscheme::impl::Extrapolation &extrapolation);

  double getValue(cplscheme::impl::Extrapolation &extrapolation, int valueID, int sampleID);
};

} // namespace testing
} // namespace precice

#pragma once

#include "precice/impl/DataContext.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the DataContext class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
class DataContextFixture {
public:
bool hasReadMapping(precice::impl::DataContext &dataContext);

bool hasWriteMapping(precice::impl::DataContext &dataContext);

void resetProvidedData(precice::impl::DataContext &dataContext);

void resetToData(precice::impl::DataContext &dataContext);

};

} // namespace testing
} // namespace precice

#pragma once

#include "precice/impl/DataContext.hpp"
#include "precice/impl/MappingContext.hpp"

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

  const impl::MappingContext mappingContext(precice::impl::DataContext &dataContext) const;

  int getFromDataID(precice::impl::DataContext &dataContext) const;

  int getToDataID(precice::impl::DataContext &dataContext) const;
};

} // namespace testing
} // namespace precice

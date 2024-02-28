#pragma once

#include <boost/range/iterator_range_core.hpp>
#include <iterator>

namespace precice::utils::ranges {

// Given a range, returns the tail of the range. Must not be empty
template <class Range>
auto tail(Range &range)
{
  return boost::make_iterator_range(++(std::begin(range)), std::end(range));
}

} // namespace precice::utils::ranges

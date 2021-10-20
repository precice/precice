#pragma once

#include "utils/span.hpp"

namespace precice {

/// Wraps a single element into a span
template <typename T>
auto refToSpan(T &element)
{
  return precice::span<T>{&element, 1};
}

} // namespace precice

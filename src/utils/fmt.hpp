/** @file
 * This file contains all required headers to consistently use fmtlib in preCICE.
 * Use this is you want to format anything in preCICE.
 */

#pragma once

#include <string>
#include <string_view>

#include "fmt/format.h"
#include "fmt/ostream.h"
#include "utils/fmtEigen.hpp"
#include "utils/fmtSTL.hpp"

namespace precice::utils {

// pass through for string only
inline std::string format_or_error(std::string_view str)
{
  return std::string{str};
}

template <class... A>
std::string format_or_error(std::string_view fmt, A &&... args)
{
  try {
    return fmt::vformat(fmt, fmt::make_format_args(args...));
  } catch (const fmt::format_error &e) {
    return std::string{"fmt_error: "} + e.what();
  }
}

} // namespace precice::utils

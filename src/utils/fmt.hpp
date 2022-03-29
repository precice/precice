/** @file
 * This file contains all required headers to consistently use fmtlib in preCICE.
 * Use this is you want to format anything in preCICE.
 */

#pragma once

#include <string>

#include "fmt/format.h"
#include "fmt/ostream.h"
#include "utils/fmtEigen.hpp"
#include "utils/fmtSTL.hpp"

template<class... A>
std::string format_or_error(A&&... args)
{
  try {
    return fmt::format(std::forward<A>(args)...);
  } catch (const fmt::format_error& e) {
    return "fmt_error: " + e.what();
  }
}

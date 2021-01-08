#pragma once

#ifdef NDEBUG

#define PRECICE_ASSERT(...) \
  {                         \
  }

#else

#include <cassert>
#include <iostream>
#include <sstream>

#include <boost/current_function.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/is_empty.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>

#include "Parallel.hpp"
#include "stacktrace.hpp"

/// Helper macro, used by assertion.
#define PRECICE_PRINT_ARGUMENT(r, out, i, elem) \
  BOOST_PP_IF(BOOST_PP_IS_EMPTY(elem), ,        \
              out << "  Argument " << i << ": " << elem << '\n';)

/// Asserts that expr evaluates to true, prints all other arguments and calls assert(false).
#define PRECICE_ASSERT(check, ...)                                                              \
  if (!(check)) {                                                                               \
    std::ostringstream oss;                                                                     \
    oss << "ASSERTION FAILED\n"                                                                 \
        << "  Location:          " << BOOST_CURRENT_FUNCTION << '\n'                            \
        << "  File:              " << __FILE__ << ":" << __LINE__ << '\n'                       \
        << "  Rank:              " << precice::utils::Parallel::getProcessRank() << '\n'        \
        << "  Failed expression: " << BOOST_PP_STRINGIZE(check)                                 \
        << '\n';                                                                                \
    BOOST_PP_SEQ_FOR_EACH_I(PRECICE_PRINT_ARGUMENT, oss, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    oss << getStacktrace() << '\n';                                                             \
    std::cerr << oss.str() << std::flush;                                                       \
    std::cout.flush();                                                                          \
    assert(false);                                                                              \
  }

#endif

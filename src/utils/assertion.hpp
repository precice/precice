#pragma once

#include <cassert>
#include <iostream>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/current_function.hpp>

#include "stacktrace.hpp"
#include "Parallel.hpp"

#ifdef NDEBUG

#define assertion(...) {}

#else

/// Helper macro, used by assertion.
#define PRINT_ARGUMENT(r, data, i, elem)                        \
  std::cerr << "  Argument " << i << ": " << elem << '\n';

/// Asserts that expr evaluates to true, prints all other arguments and calls assert(false).
#define assertion(...) if (not (BOOST_PP_VARIADIC_ELEM(0, __VA_ARGS__))) { \
    std::cerr << "ASSERTION FAILED\n" \
              << "  Location:          " << BOOST_CURRENT_FUNCTION << '\n' \
              << "  File:              " << __FILE__ << ":" << __LINE__ << '\n' \
              << "  Rank:              " << precice::utils::Parallel::getProcessRank()  << '\n' \
              << "  Failed expression: " << BOOST_PP_STRINGIZE(BOOST_PP_VARIADIC_ELEM(0, __VA_ARGS__)) \
              <<  '\n';                                            \
    BOOST_PP_SEQ_FOR_EACH_I(PRINT_ARGUMENT,, BOOST_PP_SEQ_TAIL(BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))); \
    std::cerr << getStacktrace() << '\n';                          \
    std::cerr.flush();                                                  \
    std::cout.flush();                                                  \
    assert(false);                                                      \
  }

#endif

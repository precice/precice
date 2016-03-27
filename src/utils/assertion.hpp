#pragma once

#include <cassert>
#include <iostream>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>


#ifdef NDEBUG

#define assertion(...)

#else

/// Helper macro, used by assertion.
#define PRINT_ARGUMENT(r, data, i, elem)                        \
  std::cerr << "  Argument " << i << ": " << elem << std::endl;

/// Asserts that expr evaluates to true, prints all other arguments and calls assert(false).
#define assertion(...) if (not (BOOST_PP_VARIADIC_ELEM(0, __VA_ARGS__))) { \
    std::cerr << "Assertion in " << __FILE__ << ":" << __LINE__         \
              << ", failed expression: " << BOOST_PP_STRINGIZE(BOOST_PP_VARIADIC_ELEM(0, __VA_ARGS__)) \
              <<  std::endl;                                            \
    BOOST_PP_SEQ_FOR_EACH_I(PRINT_ARGUMENT,, BOOST_PP_SEQ_TAIL(BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))); \
    std::cerr.flush();                                                  \
    std::cout.flush();                                                  \
    assert(false);                                                      \
  }

#endif

#define assertion1 assertion
#define assertion2 assertion
#define assertion3 assertion
#define assertion4 assertion
#define assertion5 assertion

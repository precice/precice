#pragma once

#include <iostream>

#ifdef NDEBUG

#define PRECICE_ASSERT(...) \
  {                         \
  }

#else

#include <cassert>
#include <sstream>

#include <boost/current_function.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/is_empty.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/variadic/size.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>

#include "fmt/format.h"
#include "utils/Parallel.hpp"
#include "utils/stacktrace.hpp"

/// Helper macro, used by TRACE
#define PRECICE_LOG_ARGUMENTS_FMT(r, data, i, elem) \
  "  " BOOST_PP_STRINGIZE(i) ": " BOOST_PP_STRINGIZE(elem) " == {}\n"

#define PRECICE_LOG_ARGUMENTS(...)                                                                                  \
  BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__),                                                                      \
              ,                                                                                                     \
              fmt::format(                                                                                          \
                  "\n" BOOST_PP_SEQ_FOR_EACH_I(PRECICE_LOG_ARGUMENTS_FMT, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)), \
                  __VA_ARGS__))

namespace precice {
namespace utils {

static constexpr char const *ASSERT_FMT =
    "ASSERTION FAILED\n"
    "Location:   {}\n"
    "File:       {}:{}\n"
    "Expression: {}\n"
    "Rank:       {}\n"
    "Arguments:  {}\n"
    "Stacktrace:\n{}\n";
}
} // namespace precice

/** Main implementation of the assertion
 * @param[in] check the expression which needs to evaluate to true for the assertion to pass
 * @param[in] args the expression which evaluates to the formatted arguments
 */
#define PRECICE_ASSERT_IMPL(check, args)                                 \
  if (!(check)) {                                                        \
    std::cerr << fmt::format(precice::utils::ASSERT_FMT,                 \
                             BOOST_CURRENT_FUNCTION, __FILE__, __LINE__, \
                             BOOST_PP_STRINGIZE(check),                  \
                             precice::utils::Parallel::getProcessRank(), \
                             args,                                       \
                             getStacktrace())                            \
              << std::flush;                                             \
    std::cout.flush();                                                   \
    assert(false);                                                       \
  }

#define PRECICE_ASSERT_IMPL_N(check, ...) \
  PRECICE_ASSERT_IMPL(check, PRECICE_LOG_ARGUMENTS(__VA_ARGS__))

#define PRECICE_ASSERT_IMPL_1(check) \
  PRECICE_ASSERT_IMPL(check, "none")

/** Asserts either a single statement or a statement followed by multiple arguments.
 *
 * This dispatches to PRECICE_ASSERT_IMPL_1 if only a check is provided, otherwise to PRECICE_ASSERT_IMPL_N
 */
#define PRECICE_ASSERT(...)                                                                                       \
  BOOST_PP_CAT(PRECICE_ASSERT_IMPL_, BOOST_PP_IF(BOOST_PP_GREATER(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__), 1), N, 1)) \
  (__VA_ARGS__)

#endif

/// Displays an error message and aborts the program independent of the build type.
/// Use to mark unreachable statements under switch or if blocks.
#define PRECICE_UNREACHABLE(message)   \
  {                                    \
    std::cerr << message << std::endl; \
    std::abort();                      \
  }

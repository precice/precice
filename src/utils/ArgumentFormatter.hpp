#pragma once

#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/vmd/is_empty.hpp>

#include "utils/fmt.hpp"

/// Helper macro, used by TRACE
#define PRECICE_LOG_ARGUMENTS_FMT(r, data, i, elem) \
  "  " BOOST_PP_STRINGIZE(i) ": " BOOST_PP_STRINGIZE(elem) " == {}\n"

#define PRECICE_LOG_ARGUMENTS(...)                                                                                  \
  BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__),                                                                      \
              "",                                                                                                   \
              ::precice::utils::format_or_error(                                                                    \
                  "\n" BOOST_PP_SEQ_FOR_EACH_I(PRECICE_LOG_ARGUMENTS_FMT, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)), \
                  __VA_ARGS__))

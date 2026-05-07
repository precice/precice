/** @file
 * This file contains all required headers to consistently use fmtlib in preCICE.
 * Use this is you want to format anything in preCICE.
 */

#pragma once

#include <format>
#include <ranges>
#include <string_view>

namespace precice::utils {

template <std::ranges::input_range R>
struct join_view {
  R const         &range;
  std::string_view sep;
  std::string_view last_sep;
};

template <std::ranges::input_range R>
join_view(R const &, std::string_view, std::string_view)
    -> join_view<R>;

template <std::ranges::input_range R>
auto join(R const &r, std::string_view sep, std::string_view last_sep)
{
  return join_view{r, sep, last_sep};
}

template <std::ranges::input_range R>
auto join(R const &r, std::string_view sep)
{
  return join_view{r, sep, sep};
}

} // namespace precice::utils

template <
    std::ranges::input_range R,
    class CharT>
struct std::formatter<precice::utils::join_view<R>, CharT> {

  std::formatter<
      std::ranges::range_value_t<R>,
      CharT>
      elem_fmt;

  constexpr auto parse(auto &ctx)
  {
    return ctx.begin();
  }

  auto format(precice::utils::join_view<R> const &jv, auto &ctx) const
  {
    auto out = ctx.out();

    auto it  = std::ranges::begin(jv.range);
    auto end = std::ranges::end(jv.range);

    if (it == end)
      return out;

    auto next = it;
    ++next;

    for (; next != end; ++it, ++next) {
      out = elem_fmt.format(*it, ctx);
      out = std::ranges::copy(jv.sep, out).out;
      ctx.advance_to(out);
    }

    out = elem_fmt.format(*it, ctx);

    if (std::ranges::size(jv.range) > 1)
      out = std::ranges::copy(jv.last_sep, out).out;

    return out;
  }
};

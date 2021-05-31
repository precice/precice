#include <array>
#include <deque>
#include <fmt/format.h>
#include <map>
#include <set>
#include <vector>

template <class T, class Allocator>
struct fmt::formatter<std::vector<T, Allocator>> : formatter<string_view> {
  /// Formats vectors like [ l, i, s, t, s ]
  template <typename FormatContext>
  auto format(const std::vector<T, Allocator> &v, FormatContext &ctx)
  {
    return format_to(ctx.out(), "[{}]", fmt::join(v, ", "));
  }
};

template <class T, class Allocator>
struct fmt::formatter<std::deque<T, Allocator>> : formatter<string_view> {
  /// Formats deques like [ l, i, s, t, s ]
  template <typename FormatContext>
  auto format(const std::deque<T, Allocator> &v, FormatContext &ctx)
  {
    return format_to(ctx.out(), "[{}]", fmt::join(v, ", "));
  }
};

template <class T, std::size_t n>
struct fmt::formatter<std::array<T, n>> : formatter<string_view> {
  /// Formats arrays like [ l, i, s, t, s ]
  template <typename FormatContext>
  auto format(const std::array<T, n> &v, FormatContext &ctx)
  {
    return format_to(ctx.out(), "[{}]", fmt::join(v, ", "));
  }
};

template <class T, class Compare, class Allocator>
struct fmt::formatter<std::set<T, Compare, Allocator>> : formatter<string_view> {
  /// Formats sets like { s, e, t, s}
  template <typename FormatContext>
  auto format(const std::set<T, Compare, Allocator> &s, FormatContext &ctx)
  {
    return format_to(ctx.out(), "{{{}}}", fmt::join(s, ", "));
  }
};

template <class T, class Compare, class Allocator>
struct fmt::formatter<std::map<T, Compare, Allocator>> : formatter<string_view> {
  /// Formats maps like tuples ( (k1, v1), (k2, v2) )
  template <typename FormatContext>
  auto format(const std::map<T, Compare, Allocator> &s, FormatContext &ctx)
  {
    return format_to(ctx.out(), "({})", fmt::join(s, ", "));
  }
};

template <class F, class S>
struct fmt::formatter<std::pair<F, S>> : formatter<string_view> {
  /// Formats pairs like tuples ( first, second )
  template <typename FormatContext>
  auto format(const std::pair<F, S> &p, FormatContext &ctx)
  {
    return format_to(ctx.out(), "({}, {})", p.first, p.second);
  }
};

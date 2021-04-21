#include <vector>
#include <fmt/format.h>

template <class T, class Allocator >
struct fmt::formatter<std::vector<T, Allocator>>: formatter<string_view> {
  /// Formats vectors like [ l, i, s, t, s ]
  template <typename FormatContext>
  auto format(const std::vector<T, Allocator>& v, FormatContext& ctx) {
    return format_to( ctx.out(), "[{}]", fmt::join(v, ", "));
  }
};

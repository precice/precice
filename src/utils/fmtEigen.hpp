#include <Eigen/Core>
#include <algorithm>
#include <fmt/format.h>
#include <sstream>
#include <utils/EigenIO.hpp>

template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
struct fmt::formatter<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>> : formatter<string_view> {
  template <typename FormatContext>
  auto format(const Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> &v, FormatContext &ctx)
  {
    std::ostringstream oss;
    oss << v.format(precice::utils::eigenio::wkt());
    const auto str = oss.str();
    std::copy(str.begin(), str.end(), ctx.out());
    return ctx.out();
  }
};

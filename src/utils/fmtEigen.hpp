#include <Eigen/Core>
#include <algorithm>
#include <format>
#include <sstream>
#include <utils/EigenIO.hpp>

// Register Eigen::Matrix as formattable type.
// We should use Eigen::DenseBase here, but this doesn't seem to work as expected.
// Maybe C++17 deduction guides will help with this?
template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
struct std::formatter<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>, char> {

  constexpr auto parse(std::format_parse_context &ctx)
  {
    return ctx.begin();
  }

  template <class FmtContext>
  FmtContext::iterator format(const Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> &v, FmtContext &ctx) const
  {
    std::ostringstream oss;
    oss << v.format(precice::utils::eigenio::wkt());
    return std::ranges::copy(oss.str(), ctx.out()).out;
  }
};

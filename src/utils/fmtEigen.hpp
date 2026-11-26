#include <Eigen/Core>
#include <algorithm>
#include <sstream>
#include "fmt/format.h"
#include "utils/EigenIO.hpp"

// Register the result of EigenDenseBase::format() as usable type
template <typename ExpressionType>
struct precice::fmt::formatter<Eigen::WithFormat<ExpressionType>> : ostream_formatter {
};

// Register Eigen::Matrix as formattable type.
// We should use Eigen::DenseBase here, but this doesn't seem to work as expected.
// Maybe C++17 deduction guides will help with this?
template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
struct precice::fmt::formatter<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>> : formatter<string_view> {
  template <typename FormatContext>
  auto format(const Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> &v, FormatContext &ctx) const
  {
    return precice::fmt::format_to(ctx.out(), "{}", v.format(precice::utils::eigenio::wkt()));
  }
};

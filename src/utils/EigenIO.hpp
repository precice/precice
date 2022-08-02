#pragma once

#include <Eigen/Core>

namespace precice {
namespace utils {
namespace eigenio {

inline Eigen::IOFormat wkt()
{
  return Eigen::IOFormat(
      Eigen::StreamPrecision,
      Eigen::DontAlignCols,
      " ", // Coeff separator
      ","  // Row separator
  );
}

inline Eigen::IOFormat debug()
{
  return wkt();
}

} // namespace eigenio
} // namespace utils
} // namespace precice

template <>
template <typename ExpressionType>
struct fmt::formatter<Eigen::WithFormat<ExpressionType>> : ostream_formatter {
};

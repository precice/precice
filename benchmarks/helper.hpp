#include <Eigen/Core>
#include <string>

#include "testing/SourceLocation.hpp"

template <int base>
Eigen::VectorXd haltonSequence(int n)
{
  Eigen::VectorXd sequence(n);
  for (int i = 0; i < n; ++i) {
    int    index    = i + 1;
    double result   = 0;
    double fraction = 1.0 / base;
    while (index > 0) {
      result += (index % base) * fraction;
      index /= base;
      fraction /= base;
    }
    sequence(i) = result;
  }
  return sequence;
}

inline Eigen::MatrixX2d generate2DHalton(int n)
{
  Eigen::MatrixX2d m(n, 2);
  m.col(0) = haltonSequence<2>(n);
  m.col(1) = haltonSequence<3>(n);
  return m;
}

inline Eigen::MatrixX3d generate3DHalton(int n)
{
  Eigen::MatrixX3d m(n, 3);
  m.col(0) = haltonSequence<2>(n);
  m.col(1) = haltonSequence<3>(n);
  m.col(2) = haltonSequence<5>(n);
  return m;
}

inline std::string benchConfig(std::string file)
{
  if (file.front() == '/') {
    file.insert(0, "/benchmarks");
  } else {
    file.insert(0, "/benchmarks/");
  }
  return precice::testing::SourceLocation + file;
}

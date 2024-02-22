#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

namespace precice::utils {

/// An accurate aggregator for doubles with usability in mind
class DoubleAggregator {
private:
  using Acc = boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::sum_kahan>>;

public:
  /// Sets the aggregator to a value
  void operator=(double d)
  {
    acc = Acc{};
    acc(d);
  }

  /// Adds a value to the aggregator
  void add(double d)
  {
    acc(d);
  }

  /// Natural version of adding a value
  void operator+=(double d)
  {
    add(d);
  }

  /// Natural version of subtracting a value
  void operator-=(double d)
  {
    add(-d);
  }

  /// Retrieves the corrected sum
  double value() const
  {
    return boost::accumulators::sum_kahan(acc);
  }

  /// Allows implicit casting to a double
  operator double() const
  {
    return value();
  }

private:
  Acc acc;
};

[[nodiscard]] inline DoubleAggregator operator+(double lhs, DoubleAggregator rhs)
{
  rhs.add(lhs);
  return rhs;
}

[[nodiscard]] inline DoubleAggregator operator+(DoubleAggregator lhs, double rhs)
{
  lhs.add(rhs);
  return lhs;
}

[[nodiscard]] inline DoubleAggregator operator-(DoubleAggregator lhs, double rhs)
{
  lhs.add(-rhs);
  return lhs;
}

} // namespace precice::utils

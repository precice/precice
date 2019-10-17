#pragma once

#include <iosfwd>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

namespace precice {
namespace utils {
namespace statistics {

/**
 * Accunulates distance measures and provides statistics based on them.
 */
class DistanceAccumulator {
public:
  /// Accumulates value
  void operator()(double value)
  {
    _acc(value);
  }

  /// Returns the minimum of all accumulated values
  double min() const
  {
    return boost::accumulators::extract::min(_acc);
  }

  /// Returns the maximum of all accumulated values
  double max() const
  {
    return boost::accumulators::extract::max(_acc);
  }

private:
  boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::min, boost::accumulators::tag::max>> _acc;
};

inline std::ostream &operator<<(std::ostream &out, const DistanceAccumulator &accumulator)
{
  return out << "min:" << accumulator.min() << " max:" << accumulator.max();
};

} // namespace statistics
} // namespace utils
} // namespace precice

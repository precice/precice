#pragma once

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <format>
#include <iosfwd>

namespace precice::utils::statistics {

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
    return empty() ? std::nan("") : boost::accumulators::extract::min(_acc);
  }

  /// Returns the maximum of all accumulated values
  double max() const
  {
    return empty() ? std::nan("") : boost::accumulators::extract::max(_acc);
  }

  /// Returns the mean of all accumulated values
  double mean() const
  {
    return empty() ? std::nan("") : boost::accumulators::extract::mean(_acc);
  }

  /// Returns the sample variance based on all accumulated values
  double variance() const
  {
    return empty() ? std::nan("") : boost::accumulators::extract::variance(_acc);
  }

  /// Returns how many values have been accumulated
  std::size_t count() const
  {
    return boost::accumulators::extract::count(_acc);
  }

  /// Returns count == 0
  bool empty() const
  {
    return count() == 0;
  }

private:
  boost::accumulators::accumulator_set<double, boost::accumulators::stats<
                                                   boost::accumulators::tag::min,
                                                   boost::accumulators::tag::max,
                                                   boost::accumulators::tag::mean,
                                                   boost::accumulators::tag::lazy_variance>>
      _acc;
};

inline std::ostream &operator<<(std::ostream &out, const DistanceAccumulator &accumulator)
{
  if (accumulator.empty()) {
    out << "empty";
  } else {
    out << "min:" << accumulator.min()
        << " max:" << accumulator.max()
        << " avg: " << accumulator.mean()
        << " var: " << accumulator.variance()
        << " cnt: " << accumulator.count();
  }
  return out;
}

} // namespace precice::utils::statistics

template <>
struct std::formatter<precice::utils::statistics::DistanceAccumulator, char> {

  constexpr auto parse(std::format_parse_context &ctx)
  {
    return ctx.begin();
  }

  template <class FmtContext>
  FmtContext::iterator format(const precice::utils::statistics::DistanceAccumulator &acc, FmtContext &ctx) const
  {
    if (acc.empty()) {
      return std::format_to(ctx.out(), "empty");
    } else {
      return std::format_to(ctx.out(), "min:{} max:{} avg:{} var:{} cnt:{}", acc.min(), acc.max(), acc.mean(), acc.variance(), acc.count());
    }
  }
};

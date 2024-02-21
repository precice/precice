#include <cmath>

namespace precice::utils {

/// A Kahan Babushka Neumaier Aggregator with usability in mind
class DoubleAggregator {
public:
  /// Sets the aggregator to a value
  void operator=(double d)
  {
    sum        = d;
    correction = 0.0;
  }

  /// Adds a value to the aggregator
  void add(double d)
  {
    double t = sum + d;
    if (std::abs(sum) >= std::abs(d)) {
      correction += (sum - t) + d;
    } else {
      correction += (d - t) + sum;
    }
    sum = t;
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

  /// Creates a new aggregator with a value added
  [[nodiscard]] DoubleAggregator operator+(double d) const
  {
    auto tmp = *this;
    tmp.add(d);
    return tmp;
  }

  /// Creates a new aggregator with a value subtracted
  [[nodiscard]] DoubleAggregator operator-(double d) const
  {
    auto tmp = *this;
    tmp.add(-d);
    return tmp;
  }

  /// Retrieves the corrected sum
  double value() const
  {
    return sum + correction;
  }

  /// Allows implicit casting to a double
  operator double() const
  {
    return value();
  }

private:
  double sum{0.0};
  double correction{0.0};
};

} // namespace precice::utils

#pragma once

#include <Eigen/Core>

namespace precice::time {

/** Result of a sampling operation which is optionally owning.
 *
 * This covers two cases:
 *
 * 1. The sampling operation doesn't need to interpolate and returns a view into an existing sample.
 * 2. The sampling operation needs to interpolate values and returns them as part of the object and exposes a view into it.
 */
class SampleResult {
public:
  /// Creates a non-owning SampleResult
  SampleResult(const Eigen::VectorXd &ref) noexcept
      : _view(ref.data(), ref.size()) {}

  /// Creates an owning SampleResult
  SampleResult(Eigen::VectorXd &&vec)
      : _owned(std::move(vec)), _view(_owned.data(), _owned.size()) {}

  // No copy
  SampleResult(const SampleResult &other)            = delete;
  SampleResult &operator=(const SampleResult &other) = delete;

  // Move is allowed
  SampleResult(SampleResult &&other)            = default;
  SampleResult &operator=(SampleResult &&other) = default;

  /// Access the values as a vector
  Eigen::Map<const Eigen::VectorXd> values() const noexcept
  {
    return _view;
  }

  /// Direct read-access to a value in the underlying view
  auto operator()(Eigen::Index index) const
  {
    return _view(index);
  }

private:
  Eigen::VectorXd                   _owned;
  Eigen::Map<const Eigen::VectorXd> _view;
};

} // namespace precice::time

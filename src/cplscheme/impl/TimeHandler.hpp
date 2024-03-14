#pragma once

#include <memory>
#include <optional>

namespace precice::cplscheme::impl {

/** Handler for stepping forward in time using time windows
 *
 * The handler respects a maximum time if it was set in the constructor.
 *
 * The goal of this class is to handle all computation around time and time steps
 * with maximizing accuracy and simplifying calling code.
 */
class TimeHandler {
public:
  TimeHandler();
  ~TimeHandler();

  /// Constructor for setting an optional maximum time
  TimeHandler(std::optional<double> maxTime);

  TimeHandler(const TimeHandler &);
  TimeHandler(TimeHandler &&);

  TimeHandler &operator=(const TimeHandler &);
  TimeHandler &operator=(TimeHandler &&);

  /// @name State altering functions
  /// @{

  /// Resets the handler to the given time
  void resetTo(double timeStart);

  /// Progress the time window by the given amount
  void progressBy(double dt);

  /// Resets the progress of the time window back to 0
  void resetProgress();

  /** Complete the time window moving the start of the time window
   *
   * The given exact time window size will be used if it is approximately equal to
   * the accumulated progress.
   */
  void completeTimeWindow(double timeWindowSize);

  /// @}
  /// @name Queries about reaching ends
  /// @{

  /** Has the end of the time window been reached?
   * This respects max time.
   */
  bool reachedEndOfWindow(double timeWindowSize) const;

  /** Has the end of the overall time been reached?
   * This is always false if no max time has been defined.
   */
  bool reachedEnd() const;

  /// @}
  /// @name Time differences
  /// @{

  /// Returns the time distance to the possibly truncated end of the current time window
  double untilWindowEnd(double timeWindowSize) const;

  /// Returns the time difference until the overall time reaches the given time
  double untilTime(double t) const;

  /** Returns the time difference until the end of the overall time.
    * This returns infinity if there is no maximum time defined.
    */
  double untilEnd() const;

  /// Returns the window progress as a double
  double windowProgress() const;

  /// @}
  /// @name Time points
  /// @{

  /// Returns the window start as a double
  double windowStart() const;

  /// Returns the current time as a double
  double time() const;

  /// @}

private:
  /// Impl to hide aggregator details
  struct Impl;

  std::unique_ptr<Impl> _impl;

  double untilProgress(double windowSize) const;
};

} // namespace precice::cplscheme::impl

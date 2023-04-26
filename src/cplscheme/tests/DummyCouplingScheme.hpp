#pragma once

#include <string>
#include <vector>
#include "../CouplingScheme.hpp"
#include "logging/Logger.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Used to test CompositionalCouplingScheme.
 */
class DummyCouplingScheme final : public CouplingScheme {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] numberIterations If 1, models an explicit coupling scheme,
   *        otherwise an implicit one.
   * @param[in] maxTimeWindows Number of time windows this DummyCouplingScheme has to perform.
   */
  DummyCouplingScheme(
      int numberIterations,
      int maxTimeWindows);

  /**
   * @brief Destructor, empty.
   */
  //virtual ~DummyCouplingScheme() {}

  /**
   * @brief
   */
  void initialize(
      double startTime,
      int    startTimesteps) override final;

  /**
   * @brief Not implemented.
   */
  void receiveResultOfFirstAdvance() override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Not implemented.
   */
  bool isInitialized() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool sendsInitializedData() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  void addComputedTime(double timeToAdd) override final
  { /* Do nothing */
  }

  /**
   * @brief
   */
  //void advance() override final;

  ChangedMeshes firstSynchronization(const ChangedMeshes &changes) override;

  void firstExchange() override;

  ChangedMeshes secondSynchronization() override;

  void secondExchange() final;

  /**
   * @brief
   */
  void finalize() override final;

  /*
   * @brief Not implemented.
   */
  std::vector<std::string> getCouplingPartners() const override final
  {
    PRECICE_ASSERT(false);
    return std::vector<std::string>();
  }

  /**
   * @brief Not implemented.
   */
  bool willDataBeExchanged(double lastSolverTimeStepSize) const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool hasDataBeenReceived() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  double getTime() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  int getTimeWindows() const override final
  {
    return _timeWindows;
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  bool hasTimeWindowSize() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  double getTimeWindowSize() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  double getNextTimeStepMaxSize() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  bool isCouplingOngoing() const override final;

  /**
   * @brief Not implemented.
   */
  bool isTimeWindowComplete() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool isActionRequired(Action action) const override final;

  bool isActionFulfilled(Action action) const override final
  {
    return true;
  }

  /**
   * @brief Not implemented.
   */
  void markActionFulfilled(Action action) override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Not implemented.
   */
  int getCheckpointTimestepInterval() const
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  void requireAction(Action action) override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Empty.
   */
  std::string printCouplingState() const override final
  {
    return std::string();
  }

  bool isImplicitCouplingScheme() const override
  {
    return _numberIterations > 1;
  }

  bool hasConverged() const override;

private:
  mutable logging::Logger _log{"cplscheme::tests::DummyCouplingScheme"};

  /// @brief Number of iterations performed per time window. 1 --> explicit.
  int _numberIterations;

  /// @brief Performed iterations in the current time window.
  int _iterations = 0;

  /// @brief Maximal number of time windows to be performed.
  int _maxTimeWindows;

  /// @brief Performed number of time windows.
  int _timeWindows = 0;

  /// @brief True, if initialize has been called.
  bool _isInitialized = false;

  /// @brief True, if timesteps are left to be performed.
  bool _isOngoing = false;

  /// @brief False, if iterations are left to be performed.
  bool _hasConverged = false;
};

} // namespace tests
} // namespace cplscheme
} // namespace precice

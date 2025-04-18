#pragma once

#include <string>
#include <vector>
#include "../CouplingScheme.hpp"
#include "logging/Logger.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::tests {

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
  // virtual ~DummyCouplingScheme() {}

  /**
   * @brief
   */
  void initialize() final override;

  void reinitialize() final override {};

  /**
   * @brief Not implemented.
   */
  bool isInitialized() const final override
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool sendsInitializedData() const final override
  {
    PRECICE_ASSERT(false);
    return false;
  }

  bool addComputedTime(double timeToAdd) final override;

  /**
   * @brief
   */
  // void advance() override final override;

  ChangedMeshes firstSynchronization(const ChangedMeshes &changes) override;

  void firstExchange() override;

  ChangedMeshes secondSynchronization() override;

  void secondExchange() final override;

  /**
   * @brief
   */
  void finalize() final override;

  /*
   * @brief Not implemented.
   */
  std::vector<std::string> getCouplingPartners() const final override
  {
    PRECICE_ASSERT(false);
    return std::vector<std::string>();
  }

  std::string localParticipant() const final override
  {
    PRECICE_ASSERT(false);
    return "unknown";
  }

  /**
   * @brief Not implemented.
   */
  bool willDataBeExchanged(double lastSolverTimeStepSize) const final override
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool hasDataBeenReceived() const final override
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  double getTime() const final override;

  double getTimeWindowStart() const final override;

  /**
   * @brief Not implemented.
   */
  int getTimeWindows() const final override
  {
    return _timeWindows;
  }

  bool hasTimeWindowSize() const final override
  {
    return true;
  }

  /**
   * @brief Not implemented.
   */
  double getTimeWindowSize() const final override
  {
    return 1.0;
  }

  /**
   * @brief Not implemented.
   */
  double getNextTimeStepMaxSize() const final override
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  bool isCouplingOngoing() const final override;

  /**
   * @brief Not implemented.
   */
  bool isTimeWindowComplete() const final override
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  bool isActionRequired(Action action) const final override;

  bool isActionFulfilled(Action action) const final override
  {
    return true;
  }

  /**
   * @brief Not implemented.
   */
  void markActionFulfilled(Action action) final override
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
  void requireAction(Action action) final override
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Empty.
   */
  std::string printCouplingState() const final override
  {
    return std::string();
  }

  bool isImplicitCouplingScheme() const override
  {
    return _numberIterations > 1;
  }

  bool hasConverged() const override;

  bool requiresSubsteps() const final override
  {
    return true;
  }

  ImplicitData implicitDataToReceive() const final override
  {
    return {};
  }

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

} // namespace precice::cplscheme::tests

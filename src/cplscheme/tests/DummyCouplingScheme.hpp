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
class DummyCouplingScheme : public CouplingScheme {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] numberIterations If 1, models and explicit coupling scheme,
   *        otherwise and implicit one.
   */
  DummyCouplingScheme(
      int numberIterations,
      int maxTimesteps);

  /**
   * @brief Destructor, empty.
   */
  virtual ~DummyCouplingScheme() {}

  /**
   * @brief
   */
  void initialize(
      double startTime,
      int    startTimesteps) override final;

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
  bool receivesInitializedData() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  void initializeData() override final
  {
    PRECICE_ASSERT(false);
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
  void advance() override final;

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
  bool willDataBeExchanged(double lastSolverTimestepLength) const override final
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
    return _timesteps;
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
  double getThisTimeWindowRemainder() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  double getNextTimestepMaxLength() const override final
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
  bool isActionRequired(const std::string &actionName) const override final;

  /**
   * @brief Not implemented.
   */
  void markActionFulfilled(const std::string &actionName) override final
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
  void requireAction(const std::string &actionName) override final
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

private:
  mutable logging::Logger _log{"cplscheme::tests::DummyCouplingScheme"};

  /// @brief Number of iterations performed per timestep. 1 --> explicit.
  int _numberIterations;

  /// @brief Performed iterations in the current timestep.
  int _iterations = 0;

  /// @brief Maximal number of timesteps to be performed.
  int _maxTimesteps;

  /// @brief Performed number of timesteps.
  int _timesteps = 0;

  /// @brief True, if initialize has been called.
  bool _isInitialized = false;

  /// @brief True, if timesteps are left to be performed.
  bool _isOngoing = false;
};

} // namespace tests
} // namespace cplscheme
} // namespace precice

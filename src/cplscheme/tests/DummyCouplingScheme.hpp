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
  virtual void initialize(
      double startTime,
      int    startTimesteps) override final;

  /**
   * @brief Not implemented.
   */
  virtual bool isInitialized() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  virtual void initializeData() override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Not implemented.
   */
  virtual void addComputedTime(double timeToAdd) override final
  { /* Do nothing */
  }

  /**
   * @brief
   */
  virtual void advance() override final;

  /**
   * @brief
   */
  virtual void finalize() override final;

  /*
   * @brief Not implemented.
   */
  virtual std::vector<std::string> getCouplingPartners() const override final
  {
    PRECICE_ASSERT(false);
    return std::vector<std::string>();
  }

  /**
   * @brief Not implemented.
   */
  virtual bool willDataBeExchanged(double lastSolverTimestepLength) const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  virtual bool hasDataBeenReceived() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  virtual double getTime() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual int getTimeWindows() const override final
  {
    return _timesteps;
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual bool hasTimeWindowSize() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  virtual double getTimeWindowSize() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual double getThisTimeWindowRemainder() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual double getNextTimestepMaxLength() const override final
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual bool isCouplingOngoing() const override final;

  /**
   * @brief Not implemented.
   */
  virtual bool isTimeWindowComplete() const override final
  {
    PRECICE_ASSERT(false);
    return false;
  }

  /**
   * @brief Not implemented.
   */
  virtual bool isActionRequired(const std::string &actionName) const override final;

  /**
   * @brief Not implemented.
   */
  virtual void markActionFulfilled(const std::string &actionName) override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Not implemented.
   */
  virtual int getCheckpointTimestepInterval() const
  {
    PRECICE_ASSERT(false);
    return 0;
  }

  /**
   * @brief Not implemented.
   */
  virtual void requireAction(const std::string &actionName) override final
  {
    PRECICE_ASSERT(false);
  }

  /**
   * @brief Empty.
   */
  virtual std::string printCouplingState() const override final
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

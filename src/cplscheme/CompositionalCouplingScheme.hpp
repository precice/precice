#pragma once

#include <list>
#include <string>
#include <vector>
#include "Constants.hpp"
#include "CouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {

/**
 * @brief Acts as one coupling scheme, but consists of several composed ones.
 *
 * The sequence of how the schemes are added by addCouplingScheme() determines
 * their execution sequence. The schemes are executed from first added to last
 * added, following the rule that implicit schemes need to be converged before
 * executing the next explicit scheme/s.
 * An example with implicit schemes I1, ..., I5 and explicit
 * schemes E1, ..., E4 is given here. The schemes are added in the sequence
 *
 *   I1, I2, E1, E2, I3, E3, I4, I5, E4
 *
 * The execution is as follows:
 *
 *   1. I1 and I2 iterate until convergence
 *   2. E1 and E2 are computed
 *   3. I3 iterates until convergence
 *   4. E3 is computed
 *   5. I4 and I5 iterate until convergence
 *   6. E4 is computed
 *
 * Steps 2 and 3 are computed directly when 1 is converged, i.e., no additional
 * solver iteration/computation is necessary. Step 3, however, will require
 * additional solver computations for every additional iteration. Similarly,
 * 4 and 5 are computed directly when 3 is converged and 6 is computed directly
 * when 5 is converged.
 *
 * The grouping into such "active" sets is repeatedly done, in order to know
 * which coupling schemes should be executed. The implicit coupling scheme
 * actions constants::actionWrite/ReadIterationCheckpoint() is used to find out
 * about whether a coupling scheme does perform iterations.
 *
 * When, in an active set of implicit coupling schemes, a scheme is converged but
 * others not, the converged scheme is put on hold until all schemes in the
 * active set are converged. This assumes that the converged state of a coupling
 * scheme is not deteriorated by the further iteration of the
 * remaining non-converged schemes.
 *
 * The execution of each scheme does not only depend on this sequence, but also
 * on how the participants are configured to be first and second in the schemes.
 * If not configured properly, a deadlock might be created.
 */
class CompositionalCouplingScheme final : public CouplingScheme {
public:
  /// Destructor, empty.
  virtual ~CompositionalCouplingScheme() {}

  /**
   * @brief Adds another coupling scheme in parallel to this scheme.
   *
   * If this coupling scheme is a normal coupling scheme, an object of
   * CompositionalCouplingScheme will be created that contains this and the
   * new scheme in parallel. If this coupling scheme is already a composed
   * scheme, the new scheme will be added as another parallel scheme.
   *
   * @return Pointer to composition of coupling schemes.
   */
  void addCouplingScheme(const PtrCouplingScheme &scheme);

  /**
   * @brief Initializes the coupling scheme and establishes a communication
   *        connection to the coupling partner.
   */
  void initialize() final override;

  void reinitialize() override final;

  /// Returns true, if any of the composed coupling schemes sendsInitializedData for this participant
  bool sendsInitializedData() const override final;

  /// Returns true, if initialize has been called.
  bool isInitialized() const final override;

  /// @copydoc cplscheme::CouplingScheme::addComputedTime()
  bool addComputedTime(double timeToAdd) final override;

  /// Exchanges data and updates the state of the coupling scheme.
  //void advance() final override;

  ChangedMeshes firstSynchronization(const ChangedMeshes &changes) override;

  void firstExchange() override;

  ChangedMeshes secondSynchronization() override;

  void secondExchange() override;

  /// Finalizes the coupling and disconnects communication.
  void finalize() final override;

  /// Returns list of all coupling partners
  std::vector<std::string> getCouplingPartners() const final override;

  /// @copydoc cplscheme::CouplingScheme::localParticipant()
  std::string localParticipant() const override final;

  /**
   * @brief Returns true, if data will be exchanged when calling advance().
   *
   * Also returns true after the last call of advance() at the end of the
   * simulation.
   *
   * @param lastSolverTimeStepSize [IN] The size of the last time step
   *        computed by the solver calling willDataBeExchanged().
   */
  bool willDataBeExchanged(double lastSolverTimeStepSize) const final override;

  /**
   * @brief checks all coupling schemes this coupling scheme is composed of.
   * @returns true, if data has been exchanged in last call of advance().
   */
  bool hasDataBeenReceived() const final override;

  /**
   * @brief Returns the currently computed time of the coupling scheme.
   *
   * This time is the minimum time of any coupling scheme in the composition.
   */
  double getTime() const final override;

  double getTimeWindowStart() const override final;

  /**
   * @brief Returns the currently computed time windows of the coupling scheme.
   *
   * The time window is the minimum time window in any coupling scheme in the composition.
   */
  int getTimeWindows() const final override;

  /**
   * @brief Returns true, if time window size is given by any of the coupling schemes in this compositional coupling scheme.
   *
   */
  bool hasTimeWindowSize() const final override;

  /**
   * @brief Returns the time window size, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid time window size is given. Check with
   * hasTimeWindowSize().
   *
   */
  double getTimeWindowSize() const final override;

  /**
   * @brief Returns the maximal size of the next time step to be computed.
   *
   * If no time step size is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   *
   * This is the minimum of all max sizes of the composed coupling schemes.
   */
  double getNextTimeStepMaxSize() const final override;

  /**
   * @brief Returns true, when the coupled simulation is still ongoing.
   *
   * As long as one composed coupling scheme is still ongoing, returns true.
   */
  bool isCouplingOngoing() const final override;

  /**
   * @brief Returns true, when the accessor can advance to the next time window.
   *
   * Only true, if all composed coupling schemes have completed the time window.
   */
  bool isTimeWindowComplete() const final override;

  /**
   * @brief Returns true, if the given action has to be performed by the accessor.
   *
   * True, if any of the composed coupling schemes requires the action.
   */
  bool isActionRequired(Action action) const final override;

  /// Returns true, if the given action has been performed by the accessor.
  bool isActionFulfilled(Action action) const final override;

  /// Tells the coupling scheme that the accessor has performed the given action.
  void markActionFulfilled(Action action) final override;

  /// Sets an action required to be performed by the accessor.
  void requireAction(Action action) final override;

  /// Returns a string representation of the current coupling state.
  std::string printCouplingState() const final override;

  /// True if one cplscheme is an implicit scheme
  bool isImplicitCouplingScheme() const final;

  /// True if the implicit scheme has converged or no implicit scheme is defined
  bool hasConverged() const final;

  /// @copydoc cplscheme::CouplingScheme::requiresSubsteps()
  bool requiresSubsteps() const override final;

  /// @copydoc cplscheme::CouplingScheme::implicitDataToReceive()
  ImplicitData implicitDataToReceive() const override final;

private:
  mutable logging::Logger _log{"cplscheme::CompositionalCouplingScheme"};

  using Schemes = std::vector<PtrCouplingScheme>;

  /// Explicit coupling schemes to be executed
  Schemes _explicitSchemes;

  /// The optional implicit scheme to be handled last
  PtrCouplingScheme _implicitScheme;

  /// Are explicit schemes on hold?
  bool _explicitOnHold = false;

  /** All schemes to run next
   *
   * This is the core of the CompositionalCouplingScheme
   *
   * All schemes start at t=0, so they all run on the first time step.
   * This is updated by finishing the complete advance calling secondExchange() using updateActiveSchemes().
   */
  std::vector<CouplingScheme *> _activeSchemes;

  /** Updates _activeSchemes to the next ones that will participate in the upcoming step in time
   * Called from secondExchange
   */
  void updateActiveSchemes();

  /// Returns all schemes in execution order, explicit as well as implicit
  std::vector<CouplingScheme *> allSchemes() const;

  /// Actions also work before initialize is called
  std::vector<CouplingScheme *> activeOrAllSchemes() const;

  /// check if time windows are compatible
  void checkCompatibleTimeWindowSizes(const CouplingScheme &impl, const CouplingScheme &expl) const;
};

} // namespace cplscheme
} // namespace precice

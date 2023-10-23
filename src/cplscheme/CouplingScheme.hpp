#pragma once

#include <map>
#include <string>
#include <vector>

#include "com/SharedPointer.hpp"
#include "precice/types.hpp"

namespace precice {
namespace cplscheme {

/**
 * @brief Interface for all coupling schemes.
 *
 * ! General description
 * A coupling scheme computes the actions to be done by the coupled participants
 * (solvers) in time. It provides interface functions to setup, advance and
 * shutdown the coupling scheme and interface functions to query the state of
 * the coupling scheme and required actions of the participants.
 *
 * ! Usage
 * -# create an object of a concrete coupling scheme class
 *    (ExplicitCouplingScheme, e.g.)
 * -# add all meshes holding data to the coupling scheme by addMesh()
 * -# configure the object by adding subclass specific information
 * -# start the coupling scheme with initialize(), where the name of the local
 *    participant, i.e. the participant using the coupling scheme object, is
 *    needed
 * -# retrieve necessary information about sent/received data and the state of
 *    the coupled simulation
 * -# query and fulfill required actions
 * -# compute data to be sent (possibly taking into account received data from
 *    initialize())
 * -# advance the coupling scheme in 4 steps:
 *    firstSynchronization, firstExchange, secondSynchronization, secondExchange;
 * -# ....
 * -# when the method isCouplingOngoing() returns false, call finalize() to
 *    stop the coupling scheme
 */
class CouplingScheme {
public:
  /// Does not define a time limit for the coupled simulation.
  static const double UNDEFINED_MAX_TIME;

  /// Does not define limit on time windows for the coupled simulation.
  static const int UNDEFINED_TIME_WINDOWS;

  /// To be used, when the time window size is determined dynamically during the coupling.
  static const double UNDEFINED_TIME_WINDOW_SIZE;

  /// To be used, when the number of max iterations is not defined (for explicit coupling).
  static const int UNDEFINED_MAX_ITERATIONS;

  /// To be used, when the number of min iterations is not defined (for explicit coupling).
  static const int UNDEFINED_MIN_ITERATIONS;

  /// Actions that are required by CouplingSchemes
  enum struct Action {
    InitializeData, ///< Is the initialization of coupling data required?
    ReadCheckpoint, ///< Is the participant required to read a previously written checkpoint?
    WriteCheckpoint ///< Is the participant required to write a checkpoint?
  };

  static std::string toString(Action action);

  CouplingScheme &operator=(CouplingScheme &&) = delete;

  virtual ~CouplingScheme() {}

  /**
   * @brief Initializes the coupling scheme and establishes a communication
   *        connection to the coupling partner. Initializes coupling data.
   *
   * @param[in] startTime starting time for coupling @BU correct?
   * @param[in] startTimeWindow counter of time window for coupling @BU correct?
   */
  virtual void initialize(
      double startTime,
      int    startTimeWindow) = 0;

  /**
   * @brief Returns whether this participant of the coupling scheme sends initialized data.
   *
   * @returns true, if this participant of the coupling scheme sends initialized data
   */
  virtual bool sendsInitializedData() const = 0;

  /// Returns true, if initialize has been called.
  virtual bool isInitialized() const = 0;

  /** @name Advancing
   *
   * Advancing the couplingscheme
   * @{
   */

  /**
   * @brief Adds newly computed time. Has to be called before every advance.
   * @param timeToAdd time to be added
   *
   * @returns true, if reaches end of the window by adding timeToAdd to time in this time step.
   */
  virtual bool addComputedTime(double timeToAdd) = 0;

  using ChangedMeshes = std::vector<MeshID>;

  /** Synchronizes mesh changes with remote participants.
   *
   * At this point, both participants may have changed the meshes.
   * Thus, we need to send local changes and receive remote changes.
   *
   * @param[in] changes MeshIDs of locally changed meshes
   *
   * @returns MeshIDs of remotely changed meshes.
   */
  virtual ChangedMeshes firstSynchronization(const ChangedMeshes &changes) = 0;

  /** Exchanges the first set of data.
   *
   * @pre \ref firstSynchronization() was called
   */
  virtual void firstExchange() = 0;

  /** Receive mesh changes from remote participants in the second step.
   *
   * At this point, the remote participant may have changed the meshes if
   * it is using a serial coupling scheme.
   * In contrast, the local participant has already communicated local changes
   * to the remote participant during @ref firstSynchronization().
   * Hence we only need to receive remote changes here.
   *
   * @note local changes are covered by \ref firstSynchronization()
   *
   * @returns MeshIDs of remotely changed meshes.
   *
   * @pre \ref firstExchange() was called
   */
  virtual ChangedMeshes secondSynchronization() = 0;

  /** Exchanges the second set of data.
   *
   * This concludes the step of the coupling scheme
   *
   * @pre \ref secondSynchronization() was called
   */
  virtual void secondExchange() = 0;

  ///@}

  /// Finalizes the coupling and disconnects communication.
  virtual void finalize() = 0;

  /// Returns list of all coupling partners.
  virtual std::vector<std::string> getCouplingPartners() const = 0;

  /**
   * @brief Returns true, if data will be exchanged when calling advance().
   *
   * Also returns true after the last call of advance() at the end of the
   * simulation.
   *
   * @param lastSolverTimeStepSize [IN] The size of the last time step
   *        computed by the solver calling willDataBeExchanged().
   */
  virtual bool willDataBeExchanged(double lastSolverTimeStepSize) const = 0;

  /// @brief Returns true, if data has been exchanged in last call of advance().
  virtual bool hasDataBeenReceived() const = 0;

  /// Returns the currently computed time of the coupling scheme.
  virtual double getTime() const = 0;

  /// Returns the currently computed time windows of the coupling scheme.
  virtual int getTimeWindows() const = 0;

  /// Returns true, if time window size is prescribed by the cpl scheme.
  virtual bool hasTimeWindowSize() const = 0;

  /**
   * @brief Returns the time window size, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid time window size is given. Check with
   * hasTimeWindowSize().
   */
  virtual double getTimeWindowSize() const = 0;

  /**
   * @brief Returns the maximal size of the next time step to be computed.
   *
   * If no time window size is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   */
  virtual double getNextTimeStepMaxSize() const = 0;

  /// Returns true, when the coupled simulation is still ongoing.
  virtual bool isCouplingOngoing() const = 0;

  /// Returns true, when the accessor can advance to the next time window.
  virtual bool isTimeWindowComplete() const = 0;

  /// Returns true, if the given action has to be performed by the accessor.
  virtual bool isActionRequired(Action action) const = 0;

  /// Returns true, if the given action has already been performed by the accessor.
  virtual bool isActionFulfilled(Action action) const = 0;

  /// Tells the coupling scheme that the accessor has performed the given action.
  virtual void markActionFulfilled(Action action) = 0;

  /// Sets an action required to be performed by the accessor.
  virtual void requireAction(Action action) = 0;

  /// Returns a string representation of the current coupling state.
  virtual std::string printCouplingState() const = 0;

  /// Returns true if the scheme or one subscheme is implicit
  virtual bool isImplicitCouplingScheme() const = 0;

  /// Returns false if the scheme is implicit and hasn't converged
  virtual bool hasConverged() const = 0;
};

} // namespace cplscheme
} // namespace precice

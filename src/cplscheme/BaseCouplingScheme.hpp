#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "Constants.hpp"
#include "CouplingData.hpp"
#include "CouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "acceleration/SharedPointer.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "impl/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {
class TXTTableWriter;
} // namespace io

namespace cplscheme {
class CouplingData;

/**
 * @brief Abstract base class for standard coupling schemes.
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
 * -# query actions and mark them as fulfilled
 * -# compute data to be sent (possibly taking into account received data from
 *    initialize())
 * -# advance the coupling scheme with advance(); where the maximum time step
 *    size (= time window size) needs to be obeyed
 * -# ....
 * -# when the method isCouplingOngoing() returns false, call finalize() to
 *    stop the coupling scheme
 */
class BaseCouplingScheme : public CouplingScheme {
public:
  enum CouplingMode { Explicit,
                      Implicit,
                      Undefined };

  BaseCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      int                           validDigits,
      std::string                   localParticipant,
      int                           maxIterations,
      CouplingMode                  cplMode,
      constants::TimesteppingMethod dtMethod);

  /**
   * @brief Getter for _sendsInitializedData
   * @returns _sendsInitializedData
   */
  bool sendsInitializedData() const override final;

  /**
   * @brief getter for _isInitialized
   * @returns true, if initialize has been called.
   */
  bool isInitialized() const override final;

  /// @copydoc cplscheme::CouplingScheme::addComputedTime()
  bool addComputedTime(double timeToAdd) override final;

  /**
   * @brief Returns true, if data will be exchanged when calling advance().
   *
   * Also returns true after the last call of advance() at the end of the
   * simulation.
   *
   * @param lastSolverTimeStepSize [IN] The size of the last time step
   *        computed by the solver calling willDataBeExchanged().
   */
  bool willDataBeExchanged(double lastSolverTimeStepSize) const override final;

  /**
   * @brief getter for _hasDataBeenReceived
   * @returns true, if data has been received in last call of advance().
   */
  bool hasDataBeenReceived() const override final;

  /**
   * @brief getter for _time
   * @returns the currently computed time of the coupling scheme.
   */
  double getTime() const override final;

  /**
   * @brief getter for _timeWindows
   * @returns the number of currently computed time windows of the coupling scheme.
   */
  int getTimeWindows() const override final;

  /**
   * @brief Function to check whether time window size is defined by coupling scheme.
   *
   * There are two reasons why a scheme might have a time window size:
   * 1) a fixed time window size is given in the scheme
   * 2) the participant received the time window size from another participant in the scheme
   *
   * @returns true, if time window size is available.
   */
  bool hasTimeWindowSize() const override final;

  /**
   * @brief Returns the time window size, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid time window size is given. Check with
   * hasTimeWindowSize().
   */
  double getTimeWindowSize() const override final;

  /**
   * @brief Returns the maximal size of the next time step to be computed.
   *
   * If no time window size is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   */
  double getNextTimeStepMaxSize() const override final;

  /// Returns true, when the coupled simulation is still ongoing.
  bool isCouplingOngoing() const override final;

  /// Returns true, when the accessor can advance to the next time window.
  bool isTimeWindowComplete() const override final;

  /// Returns true, if the given action has to be performed by the accessor.
  bool isActionRequired(Action action) const override final;

  /// Returns true, if the given action has to be performed by the accessor.
  bool isActionFulfilled(Action action) const override final;

  /// Tells the coupling scheme that the accessor has performed the given action.
  void markActionFulfilled(Action action) override final;

  /// Sets an action required to be performed by the accessor.
  void requireAction(Action action) override final;

  /**
   * @brief Returns coupling state information.
   *
   * Includes current iteration, max iterations, time, time window and action.
   */
  std::string printCouplingState() const override;

  /// Finalizes the coupling scheme.
  void finalize() override final;

  /**
   * @brief Initializes the coupling scheme.
   *
   * @param[in] startTime starting time of coupling scheme
   * @param[in] startTimeWindow starting counter of time window, from which coupling scheme starts
   */
  void initialize(double startTime, int startTimeWindow) override final;

  ChangedMeshes firstSynchronization(const ChangedMeshes &changes) override final;

  void firstExchange() override final;

  ChangedMeshes secondSynchronization() override final;

  void secondExchange() override final;

  /// Adds a measure to determine the convergence of coupling iterations.
  void addConvergenceMeasure(
      int                         dataID,
      bool                        suffices,
      bool                        strict,
      impl::PtrConvergenceMeasure measure,
      bool                        doesLogging);

  /// Set an acceleration technique.
  void setAcceleration(const acceleration::PtrAcceleration &acceleration);

  /**
   * @brief Getter for _doesFirstStep
   * @returns _doesFirstStep
   */
  bool doesFirstStep() const;

  /**
   * @returns true, if coupling scheme has any sendData
   */
  virtual bool hasAnySendData() = 0;

  /**
   * @brief Determines which data is initialized and therefore has to be exchanged during initialize.
   *
   * Calls determineInitialSend and determineInitialReceive for all send and receive data of this coupling scheme.
   */
  virtual void determineInitialDataExchange() = 0;

  /**
   * @brief Function to determine whether coupling scheme is an implicit coupling scheme
   * @returns true, if coupling scheme is implicit
   */
  bool isImplicitCouplingScheme() const override;

  /**
   * @brief Checks if the implicit cplscheme has converged
   *
   * @pre \ref doImplicitStep() or \ref receiveConvergence() has been called
   */
  bool hasConverged() const override;

protected:
  /// All send and receive data as a map "data ID -> data"
  DataMap _allData;

  /// Acceleration method to speedup iteration convergence.
  acceleration::PtrAcceleration _acceleration;

  void sendNumberOfTimeSteps(const m2n::PtrM2N &m2n, const int numberOfTimeSteps);

  void sendTimes(const m2n::PtrM2N &m2n, const Eigen::VectorXd &times);

  /**
   * @brief Sends data sendDataIDs given in mapCouplingData with communication.
   *
   * @param m2n M2N used for communication
   * @param sendData DataMap associated with sent data
   */
  void sendData(const m2n::PtrM2N &m2n, const DataMap &sendData);

  int receiveNumberOfTimeSteps(const m2n::PtrM2N &m2n);

  Eigen::VectorXd receiveTimes(const m2n::PtrM2N &m2n, int nTimeSteps);

  /**
   * @brief Receives data receiveDataIDs given in mapCouplingData with communication.
   *
   * @param m2n M2N used for communication
   * @param receiveData DataMap associated with received data
   */
  void receiveData(const m2n::PtrM2N &m2n, const DataMap &receiveData);

  /**
   * @brief Initializes storage in receiveData as zero
   *
   * @param receiveData DataMap associated with received data
   */
  void initializeWithZeroInitialData(const DataMap &receiveData);

  /**
   * @brief Adds CouplingData with given properties to this BaseCouplingScheme and returns a pointer to the CouplingData
   *
   * If CouplingData with ID of provided data already exists in coupling scheme, no duplicate is created but a pointer to the already existing CouplingData is returned.
   *
   * @param data data the CouplingData is associated with
   * @param mesh mesh the CouplingData is associated with
   * @param requiresInitialization true, if CouplingData requires initialization
   * @param exchangeSubsteps true, if CouplingData exchanges all substeps in send/recv
   * @return PtrCouplingData pointer to CouplingData owned by the CouplingScheme
   */
  PtrCouplingData addCouplingData(const mesh::PtrData &data, mesh::PtrMesh mesh, bool requiresInitialization, bool exchangeSubsteps);

  /**
   * @brief Function to determine whether coupling scheme is an explicit coupling scheme
   * @returns true, if coupling scheme is explicit
   */
  bool isExplicitCouplingScheme();

  /**
   * @brief Setter for _timeWindowSize
   * @param timeWindowSize
   */
  void setTimeWindowSize(double timeWindowSize);

  /**
   * @brief Getter for _computedTimeWindowPart
   * @returns _computedTimeWindowPart
   */
  double getComputedTimeWindowPart() const;

  /**
   * @brief Returns the time at the beginning of the current time window.
   *
   * @return time at beginning of the current time window.
   */
  double getWindowStartTime() const;

  /**
   * @brief Setter for _doesFirstStep
   */
  void setDoesFirstStep(bool doesFirstStep);

  /**
   * @brief Used to set flag after data has been received using receiveData().
   */
  void checkDataHasBeenReceived();

  /**
   * @brief Getter for _receivesInitializedData
   * @returns _receivesInitializedData
   */
  bool receivesInitializedData() const;

  /**
   * @brief Setter for _timeWindows
   *
   * Sets the computed time windows of the coupling scheme.
   * Used for testing to allow to advance in time without a coupling partner.
   *
   * @param timeWindows number of time windows
   */
  void setTimeWindows(int timeWindows);

  /**
   * @brief sends convergence to other participant via m2n
   * @param m2n used for sending
   */
  void sendConvergence(const m2n::PtrM2N &m2n);

  /**
   * @brief receives convergence from other participant via m2n
   * @param m2n used for receiving
   * @returns convergence bool
   */
  void receiveConvergence(const m2n::PtrM2N &m2n);

  /**
   * @brief perform a coupling iteration
   * @see hasConverged
   *
   * This function is called from the child classes
   */
  void doImplicitStep();

  /**
   * @brief finalizes this window's data and initializes data for next window.
   */
  void moveToNextWindow();

  /**
   * @brief used for storing all Data at end of doImplicitStep for later reference.
   */
  void storeIteration();

  /**
   * @brief Sets _sendsInitializedData, if sendData requires initialization
   * @param sendData CouplingData being checked
   */
  void determineInitialSend(DataMap &sendData);

  /**
   * @brief Sets _receivesInitializedData, if receiveData requires initialization
   * @param receiveData CouplingData being checked
   */
  void determineInitialReceive(DataMap &receiveData);

private:
  /// Coupling mode used by coupling scheme.
  CouplingMode _couplingMode = Undefined;

  mutable logging::Logger _log{"cplscheme::BaseCouplingScheme"};

  /// Maximum time being computed. End of simulation is reached, if getTime() == _maxTime
  double _maxTime;

  /// time of beginning of the current time window
  double _timeWindowStartTime = 0;

  /// Number of time windows that have to be computed. End of simulation is reached, if _timeWindows == _maxTimeWindows
  int _maxTimeWindows;

  /// number of completed time windows; _timeWindows <= _maxTimeWindows
  int _timeWindows = 0;

  /// size of time window; _timeWindowSize <= _maxTime
  double _timeWindowSize;

  /// Part of the window that is already computed; _computedTimeWindowPart <= _timeWindowSize
  double _computedTimeWindowPart = 0;

  /// Limit of iterations during one time window. Continue to next time window, if _iterations == _maxIterations.
  int _maxIterations = -1;

  /// Number of iterations in current time window. _iterations <= _maxIterations
  int _iterations = -1;

  /// Number of total iterations performed.
  int _totalIterations = -1;

  /// True, if local participant is the one starting the explicit scheme.
  bool _doesFirstStep = false;

  /// True, if _computedTimeWindowPart == _timeWindowSize and (coupling has converged or _iterations == _maxIterations)
  bool _isTimeWindowComplete = false;

  /// True, if this participant has to send initialized data.
  bool _sendsInitializedData = false;

  /// True, if this participant has to receive initialized data.
  bool _receivesInitializedData = false;

  /// True, if data has been received from other participant. Flag is used to make sure that coupling scheme is implemented and used correctly.
  bool _hasDataBeenReceived = false;

  /// True, if coupling has been initialized.
  bool _isInitialized = false;

  std::set<Action> _requiredActions;

  std::set<Action> _fulfilledActions;

  /// True if implicit scheme converged
  bool _hasConverged = false;

  /// Responsible for monitoring iteration count over time window.
  std::shared_ptr<io::TXTTableWriter> _iterationsWriter;

  /// Writes out coupling convergence within all time windows.
  std::shared_ptr<io::TXTTableWriter> _convergenceWriter;

  /// Local participant name.
  std::string _localParticipant = "unknown";

  /// Smallest number, taking validDigits into account: eps = std::pow(10.0, -1 * validDigits)
  const double _eps;

  /**
   * @brief Holds meta information to perform a convergence measurement.
   * @param data Associated data field
   * @param couplingData Coupling data history
   * @param suffices Whether this measure already suffices for convergence
   * @param strict Whether non-convergence of this measure leads to a premature end of the simulation
   * @param measure Link to the actual convergence measure
   * @param doesLogging Whether this measure is logged in the convergence file
   */
  struct ConvergenceMeasureContext {
    PtrCouplingData             couplingData;
    bool                        suffices;
    bool                        strict;
    impl::PtrConvergenceMeasure measure;
    bool                        doesLogging;

    std::string logHeader() const
    {
      return "Res" + measure->getAbbreviation() + "(" + couplingData->getDataName() + ")";
    }
  };

  /**
   * @brief All convergence measures of coupling iterations.
   *
   * Before initialization, only dataID and measure variables are filled. Then,
   * the data is fetched from send and receive data assigned to the cpl scheme.
   */
  std::vector<ConvergenceMeasureContext> _convergenceMeasures;

  /// Functions needed for initialize()

  /**
   * @brief Need to initialize receive data
   */
  virtual void initializeReceiveDataStorage() = 0;

  /**
   * @brief implements functionality for initialize in base class.
   */
  virtual void exchangeInitialData() = 0;

  /// Functions needed for advance()

  /// Exchanges the first set of data
  virtual void exchangeFirstData() = 0;

  /// Exchanges the second set of data
  virtual void exchangeSecondData() = 0;

  /**
   * @brief interface to provide accelerated data, depending on coupling scheme being used
   * @return data being accelerated
   */
  virtual const DataMap &getAccelerationData() = 0;

  /**
   * @brief If any required actions are open, an error message is issued.
   */
  void checkCompletenessRequiredActions();

  /**
   * @brief Function to check whether end of time window is reached. Does not check for convergence
   * @returns true if end time of time window is reached.
   */
  bool reachedEndOfTimeWindow();

  /**
   * @brief Initialize txt writers for iterations and convergence tracking
   */
  void initializeTXTWriters();

  /**
   * @brief Advance txt writers for iterations and convergence tracking
   */
  void advanceTXTWriters();

  /**
   * @brief Prints the coupling state
   *
   * @param timeWindows current number of completed time windows
   * @param time current time
   */
  std::string printBasicState(
      int    timeWindows,
      double time) const;

  /**
   * @brief Prints the action state
   * @returns a string representing the required actions.
   */
  std::string printActionsState() const;

  /**
   * @brief Measure whether coupling scheme has converged or not
   * @return Whether coupling scheme has converged
   */
  bool measureConvergence();

  /**
   * @brief Reset all convergence measurements after convergence
   */
  void newConvergenceMeasurements();

  /**
   * @brief Checks whether any CouplingData in dataMap requires initialization
   * @param dataMap map containing CouplingData
   * @return true, if any CouplingData in dataMap requires initialization
   */
  bool anyDataRequiresInitialization(DataMap &dataMap) const;
};
} // namespace cplscheme
} // namespace precice

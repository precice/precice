#pragma once

#include <Eigen/Core>
#include <set>
#include "Constants.hpp"
#include "CouplingData.hpp"
#include "CouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "acceleration/SharedPointer.hpp"
#include "impl/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice {
namespace cplscheme {

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
 * -# advance the coupling scheme with advance(); where the maximum timestep
 *    length (= time window size) needs to be obeyed
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
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      int                           maxIterations,
      CouplingMode                  cplMode,
      constants::TimesteppingMethod dtMethod);

  /// Adds data to be sent on data exchange and possibly be modified during coupling iterations.
  void addDataToSend(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          initialize);

  /// Adds data to be received on data exchange.
  void addDataToReceive(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          initialize);

  /// Returns true, if initialize has been called.
  virtual bool isInitialized() const
  {
    return _initializeHasBeenCalled;
  }

  /// Adds newly computed time. Has to be called before every advance.
  virtual void addComputedTime(double timeToAdd);

  /**
   * @brief Returns true, if data will be exchanged when calling advance().
   *
   * Also returns true after the last call of advance() at the end of the
   * simulation.
   *
   * @param lastSolverTimestepLength [IN] The length of the last timestep
   *        computed by the solver calling willDataBeExchanged().
   */
  virtual bool willDataBeExchanged(double lastSolverTimestepLength) const;

  /// Returns true, if data has been exchanged in last call of advance().
  virtual bool hasDataBeenExchanged() const;

  /// Returns the currently computed time of the coupling scheme.
  virtual double getTime() const;

  /// Returns the number of currently computed time windows of the coupling scheme.
  virtual int getTimeWindows() const;

  /// Returns the maximal time to be computed.
  virtual double getMaxTime() const
  {
    return _maxTime;
  }

  /// Returns the maximal number of time windows to be computed.
  virtual int getMaxTimeWindows() const
  {
    return _maxTimeWindows;
  }

  /// Returns true, if time window size is prescribed by the cpl scheme.
  virtual bool hasTimeWindowSize() const;

  /**
   * @brief Returns the time window size, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid time window size is given. Check with
   * hasTimeWindowSize().
   */
  virtual double getTimeWindowSize() const;

  /// returns list of all coupling partners
  virtual std::vector<std::string> getCouplingPartners() const;

  /**
   * @brief Returns the remaining timestep length within the current time window.
   *
   * If no time window size is prescribed by the coupling scheme, always 0.0 is
   * returned.
   */
  virtual double getThisTimeWindowRemainder() const;

  /// Returns part of the current time window that has been computed already.
  virtual double getComputedTimeWindowPart() const
  {
    return _computedTimeWindowPart;
  }

  /**
   * @brief Returns the maximal length of the next timestep to be computed.
   *
   * If no time window size is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   */
  virtual double getNextTimestepMaxLength() const; // TODO mainly used in tests. Is this function actually needed or can we drop it and only use getThisTimeWindowRemainder()?

  /// Returns true, when the coupled simulation is still ongoing.
  virtual bool isCouplingOngoing() const;

  /// Returns true, when the accessor can advance to the next time window.
  virtual bool isTimeWindowComplete() const;

  /// Returns true, if the given action has to be performed by the accessor.
  virtual bool isActionRequired(const std::string &actionName) const;

  /// Tells the coupling scheme that the accessor has performed the given action.
  virtual void markActionFulfilled(const std::string &actionName);

  /// Sets an action required to be performed by the accessor.
  virtual void requireAction(const std::string &actionName);

  /// Finalizes the coupling scheme.
  virtual void finalize();

  /**
   * @brief Initializes the coupling scheme.
   *
   * @param[in] startTime starting time of coupling scheme
   * @param[in] startTimeWindow ID of time window, from which coupling scheme starts
   */
  void initialize(double startTime, int startTimeWindow) override final;

  /**
   * @brief Initializes data with written values.
   *
   * @pre initialize() has been called.
   * @pre advance() has NOT yet been called.
   */
  void initializeData() override final;

  /**
   * @brief TODO
   */
  void advance() override final;

  /// Returns whether the solver has to evaluate the coarse or the fine model representation
  virtual bool getIsCoarseModelOptimizationActive() const
  {
    return _isCoarseModelOptimizationActive;
  }

  /// Has to be called to notify coupling scheme that coarse Model Optimization is Active
  virtual void activateCoarseModelOptimization()
  {
    _isCoarseModelOptimizationActive = true;
  }

  /**
   * @brief Sets order of predictor of interface values for first participant.
   *
   * The first participant in the implicit coupling scheme has to take some
   * initial guess for the interface values computed by the second participant.
   * In order to improve this initial guess, an extrapolation from previous
   * time windows can be performed.
   *
   * The standard predictor is of order zero, i.e., simply the converged values
   * of the last time windows are taken as initial guess for the coupling iterations.
   * Currently, an order 1 predictor is implement besides that.
   */
  void setExtrapolationOrder(int order);

  typedef std::map<int, PtrCouplingData> DataMap; // move that back to protected

  void extrapolateData(DataMap &data);

  /// Adds a measure to determine the convergence of coupling iterations.
  void addConvergenceMeasure(
      mesh::PtrData               data,
      bool                        suffices,
      int                         level,
      impl::PtrConvergenceMeasure measure);

  /// Set a coupling iteration acceleration technique.
  void setIterationAcceleration(acceleration::PtrAcceleration acceleration);

protected:
  /// Returns true, if coupling scheme is explicit
  bool isExplicitCouplingScheme()
  {
    PRECICE_CHECK(_couplingMode != Undefined, "Undefined coupling mode is not allowed.");
    return _couplingMode == Explicit;
  }

  /// Returns true, if coupling scheme is implicit
  bool isImplicitCouplingScheme()
  {
    PRECICE_CHECK(_couplingMode != Undefined, "Undefined coupling mode is not allowed.");
    return _couplingMode == Implicit;
  }

  /// Updates internal state of coupling scheme for next time window.
  void timeWindowCompleted();

  /// Receives and sets the time window size, if this participant is the one to receive
  void receiveAndSetDt();

  /// Sends the timestep length, if this participant is the one to send
  void sendDt();

  /// @return True, if local participant is the one starting the scheme.
  bool doesFirstStep() const
  {
    return _doesFirstStep;
  }

  /// Sends data sendDataIDs given in mapCouplingData with communication.
  std::vector<int> sendData(m2n::PtrM2N m2n);

  /// Receives data receiveDataIDs given in mapCouplingData with communication.
  std::vector<int> receiveData(m2n::PtrM2N m2n);

  /// Returns all data to be sent.
  DataMap &getSendData()
  {
    return _sendData;
  }

  DataMap &getReceiveData()
  {
    return _receiveData;
  }

  /// Sets the values
  CouplingData *getSendData(int dataID);

  /// Returns all data to be received with data ID as given.
  CouplingData *getReceiveData(int dataID);

  /// Sets flag to determine whether data has been exchanged in the last coupling iteration.
  void setHasDataBeenExchanged(bool hasDataBeenExchanged);

  /**
   * @brief Sets the computed time windows of the coupling scheme.
   *
   * Used for testing.
   */
  void setTimeWindows(int timeWindows)
  {
    _timeWindows = timeWindows;
  }

  /**
   * @brief Returns coupling state information.
   *
   * Includes current iteration, max iterations, time, time window and action.
   */
  virtual std::string printCouplingState() const;

  /**
   * @brief As the version without parameters, but with changed time window and time.
   *
   * This version is used by the ImplicitCouplingScheme at the moment, which
   * needs to use the last time window in the plotting when the iterations of
   * a time window are converged.
   */
  std::string printBasicState(
      int    timeWindows,
      double time) const;

  /// @return Communication device to the other coupling participant.
  m2n::PtrM2N getM2N()
  {
    PRECICE_ASSERT(_m2n);
    return _m2n;
  }

  bool hasToSendInitData()
  {
    return _hasToSendInitData;
  }

  bool hasToReceiveInitData()
  {
    return _hasToReceiveInitData;
  }

  /// Holds relevant variables to perform a convergence measurement.
  struct ConvergenceMeasure {
    mesh::PtrData               data;
    CouplingData *              couplingData;
    bool                        suffices;
    int                         level;
    impl::PtrConvergenceMeasure measure;
  };

  virtual void setupConvergenceMeasures() final;

  void newConvergenceMeasurements();

  bool measureConvergence(
      std::map<int, Eigen::VectorXd> &designSpecification);

  bool measureConvergenceCoarseModelOptimization(
      std::map<int, Eigen::VectorXd> &designSpecification);

  /**
   * @brief Sets up _dataStorage to store data values of last timestep (@BU or time window?).
   *
   * Every send data has an entry in _dataStorage. Every Entry is a vector
   * of data values with length according to the total number of values on all
   * meshes. The ordering of the data values corresponds to that in the meshes
   * and the ordering of the meshes to that in _couplingData.
   */
  void setupDataMatrices(DataMap &data);

  acceleration::PtrAcceleration getAcceleration()
  {
    return _acceleration;
  }

  int getExtrapolationOrder()
  {
    return _extrapolationOrder;
  }

  bool maxIterationsReached();

  /// TODO
  void setDeletedColumnsPPFiltering(int deletedColumnsPPFiltering)
  {
    _deletedColumnsPPFiltering = deletedColumnsPPFiltering;
  }

  /// returns number of coarse model optimization iterations in current time window
  int getIterationsCoarseOptimization(){
    return _iterationsCoarseOptimization;
  }

  /// @todo
  void initializeSendingParticipants(DataMap &dataMap);

  /// @todo
  void initializeReceivingParticipants(DataMap &dataMap);

  std::vector<ConvergenceMeasure> getConvergenceMeasures()
  {
    return _convergenceMeasures;
  }

  virtual void deactivateCoarseModelOptimization()
  {
    _isCoarseModelOptimizationActive = false;
  }

private:
  /// Communication device to the other coupling participant.
  m2n::PtrM2N _m2n;

  /// Coupling mode used by coupling scheme.
  CouplingMode _couplingMode = Undefined;

  /// Determines, if the timestep length is set by the participant.
  bool _participantSetsDt = false;

  /// Determines, if the dt length is set received from the other participant
  bool _participantReceivesDt = false;

  mutable logging::Logger _log{"cplscheme::BaseCouplingScheme"};

  /// Maximum time being computed. End of simulation is reached, if _time == _maxTime
  double _maxTime;

  /// current time; _time <= _maxTime
  double _time = 0;

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

  /// Number of coarse model optimization iterations in current time window.
  int _iterationsCoarseOptimization;

  /// Number of accumulated coarse model optimization iterations in current time window.
  int _totalIterationsCoarseOptimization = -1;

  /// TODO
  int _deletedColumnsPPFiltering = 0;

  std::vector<double> _firstResiduumNorm = {0};

  /// Extrapolation order of coupling data for first iteration of every dt.
  int _extrapolationOrder = 0;

  int _validDigits;

  /// True, if local participant is the one starting the explicit scheme.
  bool _doesFirstStep = false;

  /// True, if _computedTimeWindowPart == _timeWindowSize and (coupling has converged or _iterations == _maxIterations)
  bool _isTimeWindowComplete = false;

  /// Acceleration method to speedup iteration convergence.
  acceleration::PtrAcceleration _acceleration;

  /// To carry initData information from initialize to initData
  bool _hasToSendInitData = false;

  /// To carry initData information from initialize to initData
  bool _hasToReceiveInitData = false;

  /// True, if data has been exchanged between solvers.
  bool _hasDataBeenExchanged = false;

  /// True, if coupling has been initialized.
  bool _initializeHasBeenCalled = false;

  /// True, if initialize data has been called.
  bool _initializeDataHasBeenCalled = false;

  std::set<std::string> _actions;

  /// Map from data ID -> all send data with that ID
  DataMap _sendData;

  /// Map from data ID -> all receive data with that ID
  DataMap _receiveData;

  /// Responsible for monitoring iteration count over time window.
  std::shared_ptr<io::TXTTableWriter> _iterationsWriter;

  /// Writes out coupling convergence within all time windows.
  std::shared_ptr<io::TXTTableWriter> _convergenceWriter;

  /// First participant name.
  std::string _firstParticipant = "unknown";

  /// Second participant name.
  std::string _secondParticipant = "unknown";

  /// Local participant name.
  std::string _localParticipant = "unknown";

  /// Smallest number, taking validDigits into account: eps = std::pow(10.0, -1 * validDigits)
  const double _eps;

  /**
   * @brief All convergence measures of coupling iterations.
   *
   * Before initialization, only dataID and measure variables are filled. Then,
   * the data is fetched from send and receive data assigned to the cpl scheme.
   */
  std::vector<ConvergenceMeasure> _convergenceMeasures;

  /// Sets whether the solver evaluates the fine or the coarse model representation
  bool _isCoarseModelOptimizationActive = false;

  /// Functions needed for initialize()

  /// implements functionality needed by initialize if __couplingMode == Implicit
  virtual void initializeImplicit() = 0;

  /// implements functionality for initialize in base class.
  virtual void initializeImplementation() = 0;

  /// Functions needed for initializeData()

  /// implements functionality for initializeData in base class.
  virtual void initializeDataImpl() = 0;

  /// Functions needed for advance()

  /// implements functionality for advance in base class.
  virtual void explicitAdvance() = 0;

  /// implements functionality for advance in base class.
  virtual std::pair<bool, bool> implicitAdvance() = 0;

  /// If any required actions are open, an error message is issued.
  void checkCompletenessRequiredActions();

  /// Returns true if end time of time window is reached. Does not check for convergence
  bool subcyclingIsCompleted();

  /**
   * @brief If coupling iteration has not converged, time will be reset to beginning of
   *        window since the window has to be repeated. Iteration counters are incremented.
   *
   * @param convergence Set true, if coupling iteration in window was successful
   * @param convergenceCoarseOptimization Optional parameter, needed if manifold mapping is used
   */
  void updateTimeAndIterations(bool convergence, bool convergenceCoarseOptimization = true);

  /**
   * @brief TODO
   */
  void initializeTXTWriters();

  /**
   * @brief TODO
   */
  void advanceTXTWriters();

  /// Returns a string representing the required actions.
  std::string printActionsState() const;

  /// Implements functionality for setupConvergenceMeasures
  virtual void assignDataToConvergenceMeasure(ConvergenceMeasure* convMeasure, int dataID){
    if ((getSendData(dataID) != nullptr)) {
      convMeasure->couplingData = getSendData(dataID);
    } else {
      convMeasure->couplingData = getReceiveData(dataID);
      PRECICE_ASSERT(convMeasure->couplingData != nullptr);
    }
  }
};
} // namespace cplscheme
} // namespace precice

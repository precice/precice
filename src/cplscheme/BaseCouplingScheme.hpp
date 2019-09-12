#pragma once

#include <Eigen/Core>
#include <set>
#include "Constants.hpp"
#include "CouplingData.hpp"
#include "CouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "impl/SharedPointer.hpp"
#include "acceleration/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice
{
namespace cplscheme
{

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
 * -# query and fulfill required actions
 * -# compute data to be sent (possibly taking into account received data from
 *    initialize())
 * -# advance the coupling scheme with advance(); where the maximum timestep
 *    length needs to be obeyed
 * -# ....
 * -# when the method isCouplingOngoing() returns false, call finalize() to
 *    stop the coupling scheme
 */
class BaseCouplingScheme : public CouplingScheme
{
public:
  BaseCouplingScheme(
      double maxTime,
      int    maxTimesteps,
      double timestepLength,
      int    validDigits);

  BaseCouplingScheme(
      double                        maxTime,
      int                           maxTimesteps,
      double                        timestepLength,
      int                           validDigits,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      int                           maxIterations,
      constants::TimesteppingMethod dtMethod);

  enum CouplingMode { Explicit,
                      Implicit,
                      Undefined };

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
  //virtual PtrCouplingScheme addSchemeInParallel (PtrCouplingScheme scheme);

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
    return _isInitialized;
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

  /// Returns the currently computed timesteps of the coupling scheme.
  virtual int getTimesteps() const;

  /// Returns the maximal time to be computed.
  virtual double getMaxTime() const
  {
    return _maxTime;
  }

  /// Returns the maximal timesteps to be computed.
  virtual int getMaxTimesteps() const
  {
    return _maxTimesteps;
  }

  /// Returns true, if timestep length is prescribed by the cpl scheme.
  virtual bool hasTimestepLength() const;

  /**
   * @brief Returns the timestep length, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid timestep is given. Check with
   * hasTimestepLength().
   */
  virtual double getTimestepLength() const;

  /// returns list of all coupling partners
  virtual std::vector<std::string> getCouplingPartners() const;

  /**
   * @brief Returns the remaining timestep length of the current time step.
   *
   * If no timestep length is precribed by the coupling scheme, always 0.0 is
   * returned.
   */
  virtual double getThisTimestepRemainder() const;

  /// Returns part of the current timestep that has been computed already.
  virtual double getComputedTimestepPart() const
  {
    return _computedTimestepPart;
  }

  /**
   * @brief Returns the maximal length of the next timestep to be computed.
   *
   * If no timestep length is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   */
  virtual double getNextTimestepMaxLength() const;

  /// Returns the number of valid digits when compare times.
  int getValidDigits() const;

  /// Returns true, when the coupled simulation is still ongoing.
  virtual bool isCouplingOngoing() const;

  /// Returns true, when the accessor can advance to the next timestep.
  virtual bool isCouplingTimestepComplete() const;

  /// Returns true, if the given action has to be performed by the accessor.
  virtual bool isActionRequired(const std::string &actionName) const;

  /// Tells the coupling scheme that the accessor has performed the given action.
  virtual void performedAction(const std::string &actionName);

  /// Sets an action required to be performed by the accessor.
  virtual void requireAction(const std::string &actionName);

  /**
   * @brief Send the state of the coupling scheme to another remote scheme.
   *
   * Used in client-server approach for parallel solvers. There, the solver
   * interface does hold a coupling scheme with no data but state. The state
   * is transferred between the solver coupling scheme and the server coupling
   * scheme via sendState and receiveState.
   */
  virtual void sendState(
      com::PtrCommunication communication,
      int                   rankReceiver);

  /**
   * @brief Receive the state of the coupling scheme from another remote scheme.
   *
   * Used in client-server approach for parallel solvers. There, the solver
   * interface does hold a coupling scheme with no data but state. The state
   * is transferred between the solver coupling scheme and the server coupling
   * scheme via sendState and receiveState.
   */
  virtual void receiveState(
      com::PtrCommunication communication,
      int                   rankSender);

  /// Finalizes the coupling scheme.
  virtual void finalize();

  /// Initializes the coupling scheme.
  virtual void initialize(double startTime, int startTimestep) = 0;

  /**
   * @brief Initializes data with written values.
   *
   * @pre initialize() has been called.
   * @pre advance() has NOT yet been called.
   */
  virtual void initializeData() = 0;

  /// Returns whether the solver has to evaluate the coarse or the fine model representation
  virtual bool isCoarseModelOptimizationActive()
  {
    return _isCoarseModelOptimizationActive;
  }

  /**
   * @brief Sets order of predictor of interface values for first participant.
   *
   * The first participant in the implicit coupling scheme has to take some
   * initial guess for the interface values computed by the second participant.
   * In order to improve this initial guess, an extrapolation from previous
   * timesteps can be performed.
   *
   * The standard predictor is of order zero, i.e., simply the converged values
   * of the last timestep are taken as initial guess for the coupling iterations.
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
  /// Sets whether explicit or implicit coupling is being done.
  CouplingMode _couplingMode = Undefined;

  /// Sets whether the solver evaluates the fine or the coarse model representation
  bool _isCoarseModelOptimizationActive = false;

  /// Updates internal state of coupling scheme for next timestep.
  void timestepCompleted();

  /// Receives and set the timestep length, if this participant is the one to receive
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
  const DataMap &getSendData() const
  {
    return _sendData;
  }

  const DataMap &getReceiveData() const
  {
    return _receiveData;
  }

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

  /// Sets value for computed timestep part.
  void setComputedTimestepPart(double computedTimestepPart)
  {
    _computedTimestepPart = computedTimestepPart;
  }

  /// Sets flag to determine whether data has been exchanged in the last coupling iteration.
  void setHasDataBeenExchanged(bool hasDataBeenExchanged);

  /**
   * @brief Sets the computed time of the coupling scheme.
   *
   * Used from subclasses and when a checkpoint has been read.
   */
  void setTime(double time)
  {
    _time = time;
  }

  /**
   * @brief Sets the computed timesteps of the coupling scheme.
   *
   * Used from subclasses and when a checkpoint has been read.
   */
  void setTimesteps(int timesteps)
  {
    _timesteps = timesteps;
  }

  void setTimestepLength(double timestepLength)
  {
    _timestepLength = timestepLength;
  }

  void setIsCouplingTimestepComplete(bool isCouplingTimestepComplete)
  {
    _isCouplingTimestepComplete = isCouplingTimestepComplete;
  }

  void setIsInitialized(bool isInitialized)
  {
    _isInitialized = isInitialized;
  }

  /// If any required actions are open, an error message is issued.
  void checkCompletenessRequiredActions();

  /**
   * @brief Returns coupling state information.
   *
   * Includes current iteration, max iterations, time, timestep and action.
   */
  virtual std::string printCouplingState() const;

  /// Returns a string representing the basic state w/o actions.
  std::string printBasicState() const;

  /**
   * @brief As the version without parameters, but with changed timestep and time.
   *
   * This version is used by the ImplicitCouplingScheme at the moment, which
   * needs to use the last timestep in the plotting when the iterations of
   * a timestep are converged.
   */
  std::string printBasicState(
      int    timesteps,
      double time) const;

  /// Returns a string representing the required actions.
  std::string printActionsState() const;

  /// First participant name.
  std::string _firstParticipant = "unknown";

  /// Second participant name.
  std::string _secondParticipant = "unknown";

  /// Local participant name.
  std::string _localParticipant = "unknown";

  /// @return Communication device to the other coupling participant.
  m2n::PtrM2N getM2N()
  {
    PRECICE_ASSERT(_m2n);
    return _m2n;
  }

  void setHasToSendInitData(bool hasToSendInitData)
  {
    _hasToSendInitData = hasToSendInitData;
  }

  void setHasToReceiveInitData(bool hasToReceiveInitData)
  {
    _hasToReceiveInitData = hasToReceiveInitData;
  }

  bool hasToSendInitData()
  {
    return _hasToSendInitData;
  }

  bool hasToReceiveInitData()
  {
    return _hasToReceiveInitData;
  }

  bool participantReceivesDt()
  {
    return _participantReceivesDt;
  }

  bool participantSetsDt()
  {
    return _participantSetsDt;
  }

  /// Holds relevant variables to perform a convergence measurement.
  struct ConvergenceMeasure {
    mesh::PtrData               data;
    CouplingData *              couplingData;
    bool                        suffices;
    int                         level;
    impl::PtrConvergenceMeasure measure;
  };

  /**
   * @brief All convergence measures of coupling iterations.
   *
   * Before initialization, only dataID and measure variables are filled. Then,
   * the data is fetched from send and receive data assigned to the cpl scheme.
   */
  std::vector<ConvergenceMeasure> _convergenceMeasures;

  void setupConvergenceMeasures();

  void newConvergenceMeasurements();

  bool measureConvergence(
      std::map<int, Eigen::VectorXd> &designSpecification);

  bool measureConvergenceCoarseModelOptimization(
      std::map<int, Eigen::VectorXd> &designSpecification);

  /**
   * @brief Sets up _dataStorage to store data values of last timestep.
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

  void initializeTXTWriters();

  void advanceTXTWriters();

  void updateTimeAndIterations(bool convergence, bool convergenceCoarseOptimization = true);

  int getMaxIterations() const
  {
    return _maxIterations;
  }

  int getExtrapolationOrder()
  {
    return _extrapolationOrder;
  }

  bool maxIterationsReached();

  /// Smallest number, taking validDigists into account: eps = std::pow(10.0, -1 * validDigits)
  const double _eps;

  int _deletedColumnsPPFiltering = 0;

  /// Number of  coarse model optimization iterations in current time step.
  int _iterationsCoarseOptimization;

private:
  /// Communication device to the other coupling participant.
  m2n::PtrM2N _m2n;

  /// Determines, if the timestep length is set by the participant.
  bool _participantSetsDt = false;

  /// Determines, if the dt length is set received from the other participant
  bool _participantReceivesDt = false;

  mutable logging::Logger _log{"cplscheme::BaseCouplingScheme"};

  double _maxTime;

  int _maxTimesteps;

  /// Number of iterations in current time step.
  int _iterations = -1;

  /// Number of accumulated coarse model optimization iterations in current time step.
  int _totalIterationsCoarseOptimization = -1;

  /// Limit of iterations during one time step.
  int _maxIterations = -1;

  /// Number of total iterations performed.
  int _totalIterations = -1;

  int _timesteps = 0;

  double _timestepLength;

  double _time = 0;

  double _computedTimestepPart = 0;

  std::vector<double> _firstResiduumNorm = {0};

  /// Extrapolation order of coupling data for first iteration of every dt.
  int _extrapolationOrder = 0;

  int _validDigits;

  /// True, if local participant is the one starting the explicit scheme.
  bool _doesFirstStep = false;

  bool _isCouplingTimestepComplete = false;

  /// Acceleration method to speedup iteration convergence.
  acceleration::PtrAcceleration _acceleration;

  /// To carry initData information from initialize to initData
  bool _hasToSendInitData = false;

  /// To carry initData information from initialize to initData
  bool _hasToReceiveInitData = false;

  /// True, if data has been exchanged between solvers.
  bool _hasDataBeenExchanged = false;

  /// True, if coupling has been initialized.
  bool _isInitialized = false;

  std::set<std::string> _actions;

  /// Map from data ID -> all send data with that ID
  DataMap _sendData;

  /// Map from data ID -> all receive data with that ID
  DataMap _receiveData;

  /// Responsible for monitoring iteration count over timesteps.
  std::shared_ptr<io::TXTTableWriter> _iterationsWriter;

  /// Writes out coupling convergence within all timesteps.
  std::shared_ptr<io::TXTTableWriter> _convergenceWriter;

  int getVertexOffset(std::map<int, int> &vertexDistribution, int rank, int dim);
};
}
} // namespace precice, cplscheme

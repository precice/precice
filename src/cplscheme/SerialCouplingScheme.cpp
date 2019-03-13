#include "SerialCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "math/math.hpp"

namespace precice {
namespace cplscheme {

SerialCouplingScheme::SerialCouplingScheme
(
  double                      maxTime,
  int                         maxTimesteps,
  double                      timestepLength,
  int                         validDigits,
  const std::string&          firstParticipant,
  const std::string&          secondParticipant,
  const std::string&          localParticipant,
  m2n::PtrM2N                 m2n,
  constants::TimesteppingMethod dtMethod,
  CouplingMode                cplMode,
  int                         maxIterations)
  :
  BaseCouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits, firstParticipant,
                     secondParticipant, localParticipant, m2n, maxIterations, dtMethod)
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  assertion(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    assertion(maxIterations == 1);
  }
}

void SerialCouplingScheme::initialize
(
  double startTime,
  int    startTimestep)
{
  TRACE(startTime, startTimestep);
  assertion(not isInitialized());
  assertion(math::greaterEquals(startTime, 0.0), startTime);
  assertion(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);

  if (_couplingMode == Implicit) {
    CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
    if (not doesFirstStep()) {
      if (not _convergenceMeasures.empty()) {
        setupConvergenceMeasures(); // needs _couplingData configured
        setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
      }
      if (getPostProcessing().get() != nullptr) {
        getPostProcessing()->initialize(getSendData()); // Reserve memory, initialize
      }
    }
    else if (getPostProcessing().get() != nullptr && getPostProcessing()->getDataIDs().size()>0) {
      int dataID = *(getPostProcessing()->getDataIDs().begin());
      CHECK(getSendData(dataID) == nullptr,
            "In case of serial coupling, post-processing can be defined for "
            << "data of second participant only!");
    }
    requireAction(constants::actionWriteIterationCheckpoint());
  }

  for (DataMap::value_type & pair : getSendData()) {
    if (pair.second->initialize) {
      CHECK(not doesFirstStep(), "Only second participant can initialize data!");
      DEBUG("Initialized data to be written");
      setHasToSendInitData(true);
      break;
    }
  }

  for (DataMap::value_type & pair : getReceiveData()) {
    if (pair.second->initialize) {
      CHECK(doesFirstStep(), "Only first participant can receive initial data!");
      DEBUG("Initialized data to be received");
      setHasToReceiveInitData(true);
    }
  }

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not hasToSendInitData() && isCouplingOngoing()) {
    DEBUG("Receiving data");
    receiveAndSetDt();
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  initializeTXTWriters();
  setIsInitialized(true);
}


void SerialCouplingScheme::initializeData()
{
  TRACE();
  CHECK(isInitialized(), "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    INFO("initializeData is skipped since no data has to be initialized");
    return;
  }

  DEBUG("Initializing Data ...");

  CHECK(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
        "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  if (hasToReceiveInitData() && isCouplingOngoing() )  {
    assertion(doesFirstStep());
    DEBUG("Receiving data");
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData() && isCouplingOngoing()) {
    assertion(not doesFirstStep());
    for (DataMap::value_type & pair : getSendData()) {
      if (pair.second->oldValues.cols() == 0)
        break;
      pair.second->oldValues.col(0) = *pair.second->values;
      // For extrapolation, treat the initial value as old timestep value
      utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
    }

    // The second participant sends the initialized data to the first particpant
    // here, which receives the data on call of initialize().
    sendData(getM2N());
    receiveAndSetDt();
    // This receive replaces the receive in initialize().
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void SerialCouplingScheme::advance()
{
  TRACE(getTimesteps(), getTime());
  #ifndef NDEBUG
  for (const DataMap::value_type & pair : getReceiveData()) {
    Eigen::VectorXd& values = *pair.second->values;
    int max = values.size();
    std::ostringstream stream;
    for (int i=0; (i < max) && (i < 10); i++){
      stream << values[i] << " ";
    }
    DEBUG("Begin advance, first New Values: " << stream.str() );
  }
  #endif
  checkCompletenessRequiredActions();

  CHECK(not hasToReceiveInitData() && not hasToSendInitData(),
        "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);

  if (_couplingMode == Explicit) {
    if (math::equals(getThisTimestepRemainder(), 0.0, _eps)) {
      setIsCouplingTimestepComplete(true);
      setTimesteps(getTimesteps() + 1);
      DEBUG("Sending data...");
      sendDt();
      sendData(getM2N());

      if (isCouplingOngoing() || doesFirstStep()) {
        DEBUG("Receiving data...");
        receiveAndSetDt();
        receiveData(getM2N());
        setHasDataBeenExchanged(true);
      }
      setComputedTimestepPart(0.0);
    }
  }
  else if (_couplingMode == Implicit) {
    bool convergence = true;
    bool convergenceCoarseOptimization = true;
    bool doOnlySolverEvaluation = false;

    if (math::equals(getThisTimestepRemainder(), 0.0, _eps)) {
      DEBUG("Computed full length of iteration");
      if (doesFirstStep()) {
        sendDt();
        sendData(getM2N());
        getM2N()->receive(convergence);
        getM2N()->receive(_isCoarseModelOptimizationActive);
        if (convergence) {
          timestepCompleted();
        }
        //if (isCouplingOngoing()) {
        receiveData(getM2N());
        //}
        setHasDataBeenExchanged(true);
      }
      else {

        // get the current design specifications from the post processing (for convergence measure)
        std::map<int, Eigen::VectorXd> designSpecifications;
        if (getPostProcessing().get() != nullptr) {
          designSpecifications = getPostProcessing()->getDesignSpecification(getSendData());
        }
        // measure convergence of coupling iteration
        // measure convergence for coarse model optimization
        if(_isCoarseModelOptimizationActive){
          DEBUG("measure convergence of coarse model optimization.");
          // in case of multilevel post processing only: measure the convergence of the coarse model optimization
          convergenceCoarseOptimization = measureConvergenceCoarseModelOptimization(designSpecifications);
          // Stop, when maximal iteration count (given in config) is reached
          if (maxIterationsReached())
            convergenceCoarseOptimization = true;

          convergence = false;
          // in case of multilevel PP only: if coarse model optimization converged
          // steering the requests for evaluation of coarse and fine model, respectively
          if(convergenceCoarseOptimization){
            _isCoarseModelOptimizationActive = false;
            doOnlySolverEvaluation = true;
          }else{
            _isCoarseModelOptimizationActive = true;
          }
        }
        // measure convergence of coupling iteration
        else{
          DEBUG("measure convergence.");
          doOnlySolverEvaluation = false;

          // measure convergence of the coupling iteration,
          convergence = measureConvergence(designSpecifications);
          // Stop, when maximal iteration count (given in config) is reached
          if (maxIterationsReached())   convergence = true;
        }

        // passed by reference, modified in MM post processing. No-op for all other post-processings
        if (getPostProcessing().get() != nullptr) {
          getPostProcessing()->setCoarseModelOptimizationActive(&_isCoarseModelOptimizationActive);
        }

        // for multi-level case, i.e., manifold mapping: after convergence of coarse problem
        // we only want to evaluate the fine model for the new input, no post-processing etc..
        if (not doOnlySolverEvaluation)
        {
          // coupling iteration converged for current time step. Advance in time.
          if (convergence) {
            if (getPostProcessing().get() != nullptr) {
              _deletedColumnsPPFiltering = getPostProcessing()->getDeletedColumns();
              getPostProcessing()->iterationsConverged(getSendData());
            }
            newConvergenceMeasurements();
            timestepCompleted();

            // no convergence achieved for the coupling iteration within the current time step
          } else if (getPostProcessing().get() != nullptr) {
            getPostProcessing()->performPostProcessing(getSendData());
          }

          // extrapolate new input data for the solver evaluation in time.
          if (convergence && (getExtrapolationOrder() > 0)) {
            extrapolateData(getSendData()); // Also stores data
          }
          else { // Store data for conv. measurement, post-processing, or extrapolation
            for (DataMap::value_type& pair : getSendData()) {
              if (pair.second->oldValues.size() > 0) {
                pair.second->oldValues.col(0) = *pair.second->values;
              }
            }
            for (DataMap::value_type& pair : getReceiveData()) {
              if (pair.second->oldValues.size() > 0) {
                pair.second->oldValues.col(0) = *pair.second->values;
              }
            }
          }

          /*
          /// @todo: (Edit: Done in the solver now) need to copy coarse old values to fine old values, as first solver always sends zeros to the second solver (as pressure vals)
          //       in the serial scheme, only the sendData is registered in MM PP, we also need to register the pressure values, i.e.
          //       old fine pressure vals = old coarse pressure vals TODO: find better solution,
          //auto fineIDs = getPostProcessing()->getDataIDs();
          //for(auto id: fineIDs){
          //  std::cout<<"id: "<<id<<", fineIds.size(): "<<fineIDs.size()<<'\n';
          //  getReceiveData(id)->oldValues.column(0) = getReceiveData(id+fineIDs.size())->oldValues.column(0);
          //}
           */

        // only fine model solver evaluation is done, no PP
        } else {

          // if the coarse model problem converged within the first iteration, i.e., no post-processing at all
          // we need to register the coarse initialized data again on the fine input data,
          // otherwise the fine input data would be zero in this case, neither anything has been computed so far for the fine
          // model nor the post processing did any data registration
          // ATTENTION: assumes that coarse data is defined after fine data in same ordering.
          if(_iterationsCoarseOptimization == 1   && getPostProcessing().get() != nullptr){
            auto fineIDs = getPostProcessing()->getDataIDs();
            for (auto& fineID : fineIDs) {
              (*getSendData(fineID)->values) = getSendData(fineID+fineIDs.size()+1)->oldValues.col(0);
            }
          }
        }

        getM2N()->send(convergence);

        getM2N()->send(_isCoarseModelOptimizationActive);

        sendData(getM2N());
        
        // the second participant does not want new data in the last iteration of the last timestep
        if (isCouplingOngoing() || not convergence) {
          receiveAndSetDt();
          receiveData(getM2N());
          setHasDataBeenExchanged(true);
        }
      }

      if (not convergence) {
        DEBUG("No convergence achieved");
        requireAction(constants::actionReadIterationCheckpoint());
      }
      else {
        DEBUG("Convergence achieved");
        advanceTXTWriters();
      }
      updateTimeAndIterations(convergence, convergenceCoarseOptimization);
      setComputedTimestepPart(0.0);
    } //subcycling completed

  }
}



}}

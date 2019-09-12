#include "ParallelCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "m2n/M2N.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"
#include "math/math.hpp"

namespace precice {
namespace cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  m2n::PtrM2N           m2n,
  constants::TimesteppingMethod dtMethod,
  CouplingMode          cplMode,
  int                   maxIterations)
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
                     secondParticipant,localParticipant,m2n,maxIterations,dtMethod)
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  PRECICE_ASSERT(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    PRECICE_ASSERT(maxIterations == 1);
  }
}

void ParallelCouplingScheme::initialize
(
  double startTime,
  int    startTimestep )
{
  PRECICE_TRACE(startTime, startTimestep);
  PRECICE_ASSERT(not isInitialized());
  PRECICE_ASSERT(math::greaterEquals(startTime, 0.0), startTime);
  PRECICE_ASSERT(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);
  if (_couplingMode == Implicit) {
    PRECICE_CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
    if (not doesFirstStep()) { // second participant
      setupConvergenceMeasures(); // needs _couplingData configured
      mergeData(); // merge send and receive data for all pp calls
      setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
      if (getAcceleration().get() != nullptr) {
        getAcceleration()->initialize(getAllData()); // Reserve memory, initialize
      }
    }

    requireAction(constants::actionWriteIterationCheckpoint());
    initializeTXTWriters();
  }

  for (DataMap::value_type & pair : getSendData()) {
    if (pair.second->initialize) {
      setHasToSendInitData(true);
      break;
    }
  }
  for (DataMap::value_type & pair : getReceiveData()) {
    if (pair.second->initialize) {
      setHasToReceiveInitData(true);
      break;
    }
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);
}

void ParallelCouplingScheme::initializeData()
{
  PRECICE_TRACE("initializeData()");
  PRECICE_CHECK(isInitialized(), "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    PRECICE_INFO("initializeData is skipped since no data has to be initialized");
    return;
  }

  PRECICE_CHECK(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
        "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (hasToSendInitData()) {
      sendData(getM2N());
    }
    if (hasToReceiveInitData()) {
      receiveData(getM2N());
      setHasDataBeenExchanged(true);
    }
  }

  else { // second participant
    if (hasToReceiveInitData()) {
      receiveData(getM2N());
      setHasDataBeenExchanged(true);

      // second participant has to save values for extrapolation
      if (_couplingMode == Implicit){
        for (DataMap::value_type & pair : getReceiveData()) {
          if (pair.second->oldValues.cols() == 0)
                    break;
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
        }
      }
    }
    if (hasToSendInitData()) {
      if (_couplingMode == Implicit) {
        for (DataMap::value_type & pair : getSendData()) {
          if (pair.second->oldValues.cols() == 0)
                    break;
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
        }
      }
      sendData(getM2N());
    }
  }

  // in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void ParallelCouplingScheme::advance()
{
  if (_couplingMode == Explicit) {
    explicitAdvance();
  }
  else if (_couplingMode == Implicit) {
    implicitAdvance();
  }
}

void ParallelCouplingScheme::explicitAdvance()
{
  PRECICE_TRACE();
  checkCompletenessRequiredActions();
  PRECICE_CHECK(!hasToReceiveInitData() && !hasToSendInitData(),
        "initializeData() needs to be called before advance if data has to be initialized!");
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);

  if (math::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    setIsCouplingTimestepComplete(true);
    setTimesteps(getTimesteps() + 1);

    if (doesFirstStep()) {
      PRECICE_DEBUG("Sending data...");
      sendDt();
      sendData(getM2N());

      PRECICE_DEBUG("Receiving data...");
      receiveAndSetDt();
      receiveData(getM2N());
      setHasDataBeenExchanged(true);
    }
    else { //second participant
      PRECICE_DEBUG("Receiving data...");
      receiveAndSetDt();
      receiveData(getM2N());
      setHasDataBeenExchanged(true);

      PRECICE_DEBUG("Sending data...");
      sendDt();
      sendData(getM2N());
    }

    //both participants
    setComputedTimestepPart(0.0);
  }
}

void ParallelCouplingScheme::implicitAdvance()
{
  PRECICE_TRACE(getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  PRECICE_CHECK(!hasToReceiveInitData() && !hasToSendInitData(),
        "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  bool convergence = false;
  bool convergenceCoarseOptimization = true;
  bool doOnlySolverEvaluation = false;
  if (math::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    PRECICE_DEBUG("Computed full length of iteration");
    if (doesFirstStep()) { //First participant
      sendData(getM2N());
      getM2N()->receive(convergence);
      getM2N()->receive(_isCoarseModelOptimizationActive);
      if (convergence) {
        timestepCompleted();
      }
      receiveData(getM2N());
    }
    else { // second participant
      receiveData(getM2N());

      // get the current design specifications from the acceleration (for convergence measure)
      std::map<int, Eigen::VectorXd> designSpecifications;
      if (getAcceleration().get() != nullptr) {
        designSpecifications = getAcceleration()->getDesignSpecification(getAllData());
      }

      // measure convergence for coarse model optimization
      if(_isCoarseModelOptimizationActive){
        PRECICE_DEBUG("measure convergence of coarse model optimization.");
        // in case of multilevel acceleration only: measure the convergence of the coarse model optimization
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
        PRECICE_DEBUG("measure convergence.");
        doOnlySolverEvaluation = false;

        // measure convergence of the coupling iteration,
        convergence = measureConvergence(designSpecifications);
        // Stop, when maximal iteration count (given in config) is reached
        if (maxIterationsReached())   convergence = true;
      }

      // passed by reference, modified in MM acceleration. No-op for all other accelerations
      if (getAcceleration().get() != nullptr) {
        getAcceleration()->setCoarseModelOptimizationActive(&_isCoarseModelOptimizationActive);
      }


      // for multi-level case, i.e., manifold mapping: after convergence of coarse problem
      // we only want to evaluate the fine model for the new input, no acceleration etc..
      if (not doOnlySolverEvaluation)
      {
        if (convergence) {
          if (getAcceleration().get() != nullptr) {
            _deletedColumnsPPFiltering = getAcceleration()->getDeletedColumns();
            getAcceleration()->iterationsConverged(getAllData());
          }
          newConvergenceMeasurements();
          timestepCompleted();
        }
        else if (getAcceleration().get() != nullptr) {
          getAcceleration()->performAcceleration(getAllData());
        }

        // extrapolate new input data for the solver evaluation in time.
        if (convergence && (getExtrapolationOrder() > 0)) {
          extrapolateData(getAllData()); // Also stores data
        }
        else { // Store data for conv. measurement, acceleration, or extrapolation
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
      }else {

       // if the coarse model problem converged within the first iteration, i.e., no acceleration at all
       // we need to register the coarse initialized data again on the fine input data,
       // otherwise the fine input data would be zero in this case, neither anything has been computed so far for the fine
       // model nor the acceleration did any data registration
       // ATTENTION: assumes that coarse data is defined after fine data in same ordering.
       if (_iterationsCoarseOptimization == 1   && getAcceleration().get() != nullptr) {
         auto fineIDs = getAcceleration()->getDataIDs();
         auto& allData = getAllData();
         for(auto& fineID : fineIDs) {
           *allData.at( fineID )->values = allData.at( fineID+fineIDs.size() )->oldValues.col(0);
         }
       }
     }

      getM2N()->send(convergence);
      getM2N()->send(_isCoarseModelOptimizationActive);

      sendData(getM2N());
    }

    // both participants
    if (not convergence) {
      PRECICE_DEBUG("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
    }
    else {
      PRECICE_DEBUG("Convergence achieved");
      advanceTXTWriters();
    }
    updateTimeAndIterations(convergence, convergenceCoarseOptimization);
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  } // subcycling complete
}



void ParallelCouplingScheme::mergeData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration." );
  PRECICE_ASSERT(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(), getSendData().end());
  _allData.insert(getReceiveData().begin(), getReceiveData().end());
}

}}

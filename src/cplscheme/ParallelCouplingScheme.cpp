#include "ParallelCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

logging::Logger ParallelCouplingScheme::_log("precice::cplscheme::ParallelCouplingScheme" );

ParallelCouplingScheme::ParallelCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  m2n::M2N::SharedPointer           m2n,
  constants::TimesteppingMethod dtMethod,
  CouplingMode          cplMode,
  int                   maxIterations)
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
                     secondParticipant,localParticipant,m2n,maxIterations,dtMethod),
  _allData ()
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  assertion(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    assertion(maxIterations == 1);
  }
}

void ParallelCouplingScheme::initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);
  if (_couplingMode == Implicit) {
    preciceCheck(not getSendData().empty(), "initialize()", "No send data configured! Use explicit scheme for one-way coupling.");
    if (not doesFirstStep()) { // second participant
      setupConvergenceMeasures(); // needs _couplingData configured
      mergeData(); // merge send and receive data for all pp calls
      setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
      if (getPostProcessing().get() != nullptr) {
        preciceCheck(getPostProcessing()->getDataIDs().size()==2 ,"initialize()",
                     "For parallel coupling, the number of coupling data vectors has to be 2, not: "
                     << getPostProcessing()->getDataIDs().size());
        getPostProcessing()->initialize(getAllData()); // Reserve memory, initialize
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
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
               "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceCheck(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
               "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (hasToSendInitData()) {
      getM2N()->startSendPackage(0);
      sendData(getM2N());
      getM2N()->finishSendPackage();
    }
    if (hasToReceiveInitData()) {
      getM2N()->startReceivePackage(0);
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
  }

  else { // second participant
    if (hasToReceiveInitData()) {
      getM2N()->startReceivePackage(0);
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
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
      getM2N()->startSendPackage(0);
      sendData(getM2N());
      getM2N()->finishSendPackage();
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
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);

  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    setIsCouplingTimestepComplete(true);
    setTimesteps(getTimesteps() + 1);

    if (doesFirstStep()) {
      DEBUG("Sending data...");
      getM2N()->startSendPackage(0);
      sendDt();
      sendData(getM2N());
      getM2N()->finishSendPackage();

      DEBUG("Receiving data...");
      getM2N()->startReceivePackage(0);
      receiveAndSetDt();
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
    else { //second participant
      DEBUG("Receiving data...");
      getM2N()->startReceivePackage(0);
      receiveAndSetDt();
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      setHasDataBeenExchanged(true);

      DEBUG("Sending data...");
      getM2N()->startSendPackage(0);
      sendDt();
      sendData(getM2N());
      getM2N()->finishSendPackage();
    }

    //both participants
    setComputedTimestepPart(0.0);
  }
}

void ParallelCouplingScheme::implicitAdvance()
{
  preciceTrace("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  bool convergence = false;
  bool convergenceCoarseOptimization = true;
  bool doOnlySolverEvaluation = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    DEBUG("Computed full length of iteration");
    if (doesFirstStep()) { //First participant
      getM2N()->startSendPackage(0);
      sendData(getM2N());
      getM2N()->finishSendPackage();
      getM2N()->startReceivePackage(0);
      getM2N()->receive(convergence);
      getM2N()->startReceivePackage(0);
      getM2N()->receive(convergenceCoarseOptimization);
      getM2N()->startReceivePackage(0);
      getM2N()->receive(_isCoarseModelOptimizationActive);
      if (convergence) {
        timestepCompleted();
      }
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
    }
    else { // second participant

      // DIRTY HACK for manifold mapping
      auto fine_IDs = getPostProcessing()->getDataIDs();
      auto& allData = getAllData();
      // TODO: make this error safe ... should be receive data or something.... and then check whether it is fine data
      Eigen::VectorXd storeOldInput = (*allData.at( fine_IDs[1] )->values);

      getM2N()->startReceivePackage(0);
      receiveData(getM2N());
      getM2N()->finishReceivePackage();

      // get the current design specifications from the post processing (for convergence measure)
      std::map<int, Eigen::VectorXd> designSpecifications;
      if (getPostProcessing().get() != nullptr) {
        designSpecifications = getPostProcessing()->getDesignSpecification(getAllData());
      }

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

          // write back the correct input for the fine model as registered in the MM, as this was overwritten by receive
          // the fine model evaluation needs the same input as the coarse model, just before convegrence
          // precice's send and receive of data automatically overrides the data, as there is different data for the respective
          // id present at the first participants side.
          // TODO: make this error safe ... should be receive data or something.... and then check whether it is fine data
          (*allData.at( fine_IDs[1] )->values) = storeOldInput;

          // reset the convergence measures for the coarse model optimization
          newConvergenceMeasurements(1);
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
        if (convergence) {
          if (getPostProcessing().get() != nullptr) {
            _deletedColumnsPPFiltering = getPostProcessing()->getDeletedColumns();
            getPostProcessing()->iterationsConverged(getAllData());
          }
          newConvergenceMeasurements();
          timestepCompleted();
        }
        else if (getPostProcessing().get() != nullptr) {
          getPostProcessing()->performPostProcessing(getAllData());
        }

        // extrapolate new input data for the solver evaluation in time.
        if (convergence && (getExtrapolationOrder() > 0)) {
          extrapolateData(getAllData()); // Also stores data
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
      }else {

       // if the coarse model problem converged within the first iteration, i.e., no post-processing at all
       // we need to register the coarse initialized data again on the fine input data,
       // otherwise the fine input data would be zero in this case, neither anything has been computed so far for the fine
       // model nor the post processing did any data registration
       // ATTENTION: assumes that coarse data is defined after fine data in same ordering.
       if (_iterationsCoarseOptSum == 1   && getPostProcessing().get() != nullptr) {
         auto fineIDs = getPostProcessing()->getDataIDs();
         auto& allData = getAllData();
         for(auto& fineID : fineIDs) {
           *allData.at( fineID )->values = allData.at( fineID+fineIDs.size() )->oldValues.col(0);
         }
       }
     }

      // send convergence
      getM2N()->startSendPackage(0);
      getM2N()->send(convergence);
      getM2N()->startSendPackage(0);
      getM2N()->send(convergenceCoarseOptimization);
      // send flag that indicates which model is evaluated
      getM2N()->startSendPackage(0);
      getM2N()->send(_isCoarseModelOptimizationActive);

      sendData(getM2N());
      getM2N()->finishSendPackage();
    }

    // both participants
    if (not convergence) {
      DEBUG("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
    }
    else {
      DEBUG("Convergence achieved");
      advanceTXTWriters();
    }
    updateTimeAndIterations(convergence, convergenceCoarseOptimization);
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  } // subcycling complete
}



void ParallelCouplingScheme::mergeData()
{
  preciceTrace("mergeData()");
  assertion(!doesFirstStep(), "Only the second participant should do the post processing." );
  assertion(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(), getSendData().end());
  _allData.insert(getReceiveData().begin(), getReceiveData().end());
}

}}

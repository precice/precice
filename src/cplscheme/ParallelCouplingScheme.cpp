#include "ParallelCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log ParallelCouplingScheme::_log("precice::cplscheme::ParallelCouplingScheme" );

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
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
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
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
    if (hasToSendInitData()) {
      if (_couplingMode == Implicit) {
        for (DataMap::value_type & pair : getSendData()) {
          if (pair.second->oldValues.cols() == 0)
                    break;
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
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
      preciceDebug("Sending data...");
      getM2N()->startSendPackage(0);
      sendDt();
      sendData(getM2N());
      getM2N()->finishSendPackage();

      preciceDebug("Receiving data...");
      getM2N()->startReceivePackage(0);
      receiveAndSetDt();
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
    else { //second participant
      preciceDebug("Receiving data...");
      getM2N()->startReceivePackage(0);
      receiveAndSetDt();
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      setHasDataBeenExchanged(true);

      preciceDebug("Sending data...");
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
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  bool convergence = false;
  bool convergenceCoarseOptimization = true;
  bool doOnlySolverEvaluation = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()) { //First participant
      std::cout<<" +++ First Participant, start sendData() (sends pressure values) +++"<<std::endl;
      getM2N()->startSendPackage(0);
      sendData(getM2N());
      getM2N()->finishSendPackage();
      std::cout<<" +++ First Participant, finish sendData() +++"<<std::endl;
      getM2N()->startReceivePackage(0);
      getM2N()->receive(convergence);
      getM2N()->startReceivePackage(0);
      getM2N()->receive(_isCoarseModelOptimizationActive);
      if (convergence) {
        timestepCompleted();
      }
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      std::cout<<" +++ First Participant, finish receiveData() (receives displ values) +++"<<std::endl;
    }
    else { // second participant
      std::cout<<" +++ Second Participant, start receiveData() (receives pressure values)  +++"<<std::endl;
      getM2N()->startReceivePackage(0);
      receiveData(getM2N());
      getM2N()->finishReceivePackage();
      std::cout<<" +++ Second Participant, finish receiveData() +++"<<std::endl;

      // -------- NEW
      std::cout<<"\n ### coarse model opt active (1): "<<_isCoarseModelOptimizationActive<<std::endl;

      auto designSpecifications = getPostProcessing()->getDesignSpecification(getAllData());

      // measure convergence for coarse model optimization
      if(_isCoarseModelOptimizationActive){

        /*
        for(auto elem : _allData){
          std::cout<<"\n data with ID: "<<elem.first<<"\n"<<(*elem.second->values)<<std::endl;
        }
        */

        preciceDebug("measure convergence of coarse model optimization.");
        // in case of multilevel post processing only: measure the convergence of the coarse model optimization
        convergenceCoarseOptimization = measureConvergenceCoarseModelOptimization(designSpecifications);
        // Stop, when maximal iteration count (given in config) is reached
        if (maxIterationsReached())   convergenceCoarseOptimization = true;

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
      //if(not _isCoarseModelOptimizationActive && convergenceCoarseOptimization){
      else{
        preciceDebug("measure convergence.");
        doOnlySolverEvaluation = false;

        /*
        for(auto elem : _allData){
          std::cout<<"\n data with ID: "<<elem.first<<"\n"<<(*elem.second->values)<<std::endl;
        }
        */

        // measure convergence of the coupling iteration,
        convergence = measureConvergence(designSpecifications);
        // Stop, when maximal iteration count (given in config) is reached
        if (maxIterationsReached())   convergence = true;
      }

      std::cout<<" ### coarse model opt active (2): "<<_isCoarseModelOptimizationActive<<std::endl;
      std::cout<<" ### coarse model converged: "<<convergenceCoarseOptimization<<std::endl;
      std::cout<<" ### fine model converged: "<<convergence<<"\n"<<std::endl;

      // passed by reference, modified in MM post processing. No-op for all other post-processings
      getPostProcessing()->setCoarseModelOptimizationActive(&_isCoarseModelOptimizationActive);

      // -------- end NEW

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

        /**
        std::cout << "before sending, after PP" << std::endl;
        for (auto elem : _allData) {
          std::cout << "\n data with ID: " << elem.first << "\n" << (*elem.second->values) << std::endl;
        }
        */

        //if (isCouplingOngoing()) {
        if (convergence && (getExtrapolationOrder() > 0)) {
          extrapolateData(getAllData()); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          for (DataMap::value_type& pair : getSendData()) {
            if (pair.second->oldValues.size() > 0) {
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
          for (DataMap::value_type& pair : getReceiveData()) {
            if (pair.second->oldValues.size() > 0) {
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
        }
      }

      getM2N()->startSendPackage(0);
      getM2N()->send(convergence);

      getM2N()->startSendPackage(0);
      getM2N()->send(_isCoarseModelOptimizationActive);

      sendData(getM2N());
      getM2N()->finishSendPackage();
      std::cout<<" +++ Second Participant, finish sendData() (send displ values)  +++"<<std::endl;
    }

    // both participants
    if (not convergence) {
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
    }
    else {
      preciceDebug("Convergence achieved");
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
  assertion1(!doesFirstStep(), "Only the second participant should do the post processing." );
  assertion1(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(), getSendData().end());
  _allData.insert(getReceiveData().begin(), getReceiveData().end());
}

}}

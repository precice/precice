#include "ParallelCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"

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
  com::PtrCommunication communication,
  int                   maxIterations,
  constants::TimesteppingMethod dtMethod )
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
			 secondParticipant,localParticipant,communication,maxIterations,dtMethod),
  _allData ()
{}



// identisch, bis auf das if(not doesFirstStep()). 
void ParallelCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  assertion(getCommunication()->isConnected());
  preciceCheck(not getSendData().empty(), "initialize()",
               "No send data configured!");
  setTime(startTime);
  setTimesteps(startTimestep);
  if (not doesFirstStep()) { // second participant
    // NOTE: Block sollte nur ausgeführt werden für implizites Schema
    setupConvergenceMeasures(); // needs _couplingData configured
    mergeData(); // merge send and receive data for all pp calls
    setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
    if (getPostProcessing().get() != NULL) {
      preciceCheck(getPostProcessing()->getDataIDs().size()==2 ,"initialize()",
		   "For parallel coupling, the number of coupling data vectors has to be 2, not: "
		   << getPostProcessing()->getDataIDs().size());
      getPostProcessing()->initialize(getAllData()); // Reserve memory, initialize
    }
  }

  requireAction(constants::actionWriteIterationCheckpoint());
  initializeTXTWriters();

  foreach (DataMap::value_type & pair, getSendData()) {
    if (pair.second->initialize) {
      setHasToSendInitData(true);
      break;
    }
  }
  foreach (DataMap::value_type & pair, getReceiveData()) {
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


// identisch für explizit und implizit
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
    if(hasToSendInitData()) {
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
    if (hasToReceiveInitData()) {
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
  }

  else { // second participant
    if (hasToReceiveInitData()) {
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);

      // second participant has to save values for extrapolation
      if (getExtrapolationOrder() > 0) {
        foreach (DataMap::value_type & pair, getReceiveData()){
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
    if (hasToSendInitData()) {
      if (getExtrapolationOrder() > 0) {
        foreach (DataMap::value_type & pair, getSendData()) {
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
  }

  // in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
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

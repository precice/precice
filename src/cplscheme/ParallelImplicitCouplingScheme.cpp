// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelImplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "Constants.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"

namespace precice {
namespace cplscheme {

tarch::logging::Log ParallelImplicitCouplingScheme::
_log("precice::cplscheme::ParallelImplicitCouplingScheme" );

ParallelImplicitCouplingScheme:: ParallelImplicitCouplingScheme
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
  ParallelCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
			 secondParticipant,localParticipant,communication,maxIterations,dtMethod)
{
  couplingMode = Implicit;
}

void ParallelImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
	       "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)) {
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()) { //First participant
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
      getCommunication()->startReceivePackage(0);
      getCommunication()->receive(convergence, 0);
      if (convergence) {
        timestepCompleted();
      }
      if (isCouplingOngoing()) {
        receiveData(getCommunication());
      }
      getCommunication()->finishReceivePackage();
    }
    else { // second participant

      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();

      convergence = measureConvergence();

      assertion2((getIterations() <= getMaxIterations()) || (getMaxIterations() == -1),
		 getIterations(), getMaxIterations());
      // Stop, when maximal iteration count (given in config) is reached
      if (getIterations() == getMaxIterations()-1) {
        convergence = true;
      }
      if (convergence) {
        if (getPostProcessing().get() != NULL) {
          getPostProcessing()->iterationsConverged(getAllData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (getPostProcessing().get() != NULL) {
        getPostProcessing()->performPostProcessing(getAllData());
      }
      getCommunication()->startSendPackage(0);
      getCommunication()->send(convergence, 0);

      if (isCouplingOngoing()) {
        if (convergence && (getExtrapolationOrder() > 0)){
          extrapolateData(getAllData()); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()) {
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()) {
            if (pair.second->oldValues.size() > 0) {
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
        }
        sendData(getCommunication());
      }
      getCommunication()->finishSendPackage();
    }

    // both participants
    if (not convergence) {
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
      increaseIterations();
      increaseTotalIterations();
      // The computed timestep part equals the timestep length, since the
      // timestep remainder is zero. Subtract the timestep length do another
      // coupling iteration.
      assertion(tarch::la::greater(getComputedTimestepPart(), 0.0));
      setTime(getTime() - getComputedTimestepPart());
    }
    else {
      preciceDebug("Convergence achieved");
      getIterationsWriter().writeData("Timesteps", getTimesteps());
      getIterationsWriter().writeData("Total Iterations", getTotalIterations());
      getIterationsWriter().writeData("Iterations", getIterations());
      int converged = getIterations() < getMaxIterations() ? 1 : 0;
      getIterationsWriter().writeData("Convergence", converged);
      setIterations(0);
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  } // subcycling complete

  // When the iterations of one timestep are converged, the old time, timesteps,
  // and iteration should be plotted, and not the 0th of the new timestep. Thus,
  // the plot values are only updated when no convergence was achieved.
  if (not convergence) {
    setTimestepToPlot(getTimesteps());
    setTimeToPlot(getTime());
    setIterationToPlot(getIterations());
  }
  else {
    increaseIterationToPlot();
  }
}

std::string ParallelImplicitCouplingScheme:: printCouplingState() const
{
  std::ostringstream os;
  os << "it " << _iterationToPlot; //_iterations;
  if(getMaxIterations() != -1 ){
    os << " of " << getMaxIterations();
  }
  os << " | " << printBasicState(_timestepToPlot, _timeToPlot) << " | " << printActionsState();
  return os.str();
}

}} // namespace precice, cplscheme

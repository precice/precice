// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialImplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "Constants.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log SerialImplicitCouplingScheme::
_log("precice::cplscheme::SerialImplicitCouplingScheme" );

SerialImplicitCouplingScheme:: SerialImplicitCouplingScheme
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
  ImplicitCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
			 secondParticipant,localParticipant,communication,maxIterations,dtMethod)
{}


void SerialImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  //tmp debug -> here wrong values
  foreach (DataMap::value_type & pair, getReceiveData()){
    utils::DynVector& values = *pair.second->values;
    preciceDebug("Begin advance, New Values: " << values);
  }
  checkCompletenessRequiredActions();

  preciceCheck(not hasToReceiveInitData() && not hasToSendInitData(), "advance()",
	       "initializeData() needs to be called before advance if data has to be initialized!");
  
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()){
      getCommunication()->startSendPackage(0);
      if (participantSetsDt()){
        preciceDebug("sending timestep length of " << getComputedTimestepPart());
        getCommunication()->send(getComputedTimestepPart(), 0);
      }
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
      getCommunication()->startReceivePackage(0);
      getCommunication()->receive(convergence, 0);
      if (convergence){
        timestepCompleted();
      }
      if (isCouplingOngoing()){
        receiveData(getCommunication());
      }
      getCommunication()->finishReceivePackage();
    }
    else {
      convergence = measureConvergence();
      assertion2((getIterations() <= getMaxIterations()) || (getMaxIterations() == -1),
                 getIterations(), getMaxIterations());
      // Stop, when maximal iteration count (given in config) is reached
      if (getIterations() == getMaxIterations()-1){
        convergence = true;
      }
      if (convergence){
        if (getPostProcessing().get() != NULL){
          getPostProcessing()->iterationsConverged(getSendData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (getPostProcessing().get() != NULL){
        getPostProcessing()->performPostProcessing(getSendData());
      }
      getCommunication()->startSendPackage(0);
      getCommunication()->send(convergence, 0);
      if (isCouplingOngoing()){
        if (convergence && (getExtrapolationOrder() > 0)){
          extrapolateData(getSendData()); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()){
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()){
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
        }
        sendData(getCommunication());
        getCommunication()->finishSendPackage();
        getCommunication()->startReceivePackage(0);
        if (participantReceivesDt()){
          double dt = UNDEFINED_TIMESTEP_LENGTH;
          getCommunication()->receive(dt, 0);
          assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
          setTimestepLength(dt);
        }
        receiveData(getCommunication());
        getCommunication()->finishReceivePackage();
      }
      else {
        getCommunication()->finishSendPackage();
      }
    }
    
    if (not convergence){
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
  }//subcycling completed

  // When the iterations of one timestep are converged, the old time, timesteps,
  // and iteration should be plotted, and not the 0th of the new timestep. Thus,
  // the plot values are only updated when no convergence was achieved.
  if (not convergence){
    setTimestepToPlot(getTimesteps());
    setTimeToPlot(getTime());
    setIterationToPlot(getIterations());
  }
  else {
    increaseIterationToPlot();
  }
  
}

}} // namespace precice, cplscheme

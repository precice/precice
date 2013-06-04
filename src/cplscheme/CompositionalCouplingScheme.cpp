// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CompositionalCouplingScheme.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log CompositionalCouplingScheme::
   _log("precice::cplscheme::CompositionalCouplingScheme");

CompositionalCouplingScheme:: CompositionalCouplingScheme
(
  PtrCouplingScheme initialCouplingScheme )
:
  CouplingScheme(
    initialCouplingScheme->CouplingScheme::getMaxTime(),
    initialCouplingScheme->getMaxTimesteps(),
    initialCouplingScheme->getTimestepLength(),
    initialCouplingScheme->getValidDigits())
{}

void CompositionalCouplingScheme:: addParallelCouplingScheme
(
  PtrCouplingScheme couplingScheme )
{
  preciceTrace("addParallelCouplingScheme()");
  _couplingSchemes.push_back(couplingScheme);
  assertion2(couplingScheme->getMaxTimesteps() == getMaxTimesteps(),
             couplingScheme->getMaxTimesteps(), getMaxTimesteps());
  assertion2(couplingScheme->getTimestepLength() == getTimestepLength(),
             couplingScheme->getTimestepLength(), getTimestepLength());
  assertion2(couplingScheme->getValidDigits() == getValidDigits(),
             couplingScheme->getValidDigits(), getValidDigits());
}

void CompositionalCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2("initialize()", startTime, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->initialize(startTime, startTimestep);
    setHasDataBeenExchanged(couplingScheme->hasDataBeenExchanged() || hasDataBeenExchanged());
    if (getTimestepLength() > couplingScheme->getTimestepLength()){
      setTimestepLength(couplingScheme->getTimestepLength());
    }
  }
  setIsInitialized(true);
}

void CompositionalCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->initializeData();
    setHasDataBeenExchanged(couplingScheme->hasDataBeenExchanged() || hasDataBeenExchanged());
  }
}


void CompositionalCouplingScheme:: addComputedTime
(
  double timeToAdd )
{
  preciceTrace1("addComputedTime()", timeToAdd);
  CouplingScheme::addComputedTime(timeToAdd);
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->addComputedTime(timeToAdd);
    assertion2(couplingScheme->getTime() == getTime(),
               couplingScheme->getTime(), getTime());
  }
}

void CompositionalCouplingScheme:: advance()
{
  preciceTrace("advance()");
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(true); // Merged with coupling schemes' states
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->advance();
    setHasDataBeenExchanged(
      couplingScheme->hasDataBeenExchanged() || hasDataBeenExchanged());
    setIsCouplingTimestepComplete(
      couplingScheme->isCouplingTimestepComplete() && isCouplingTimestepComplete());
    if (getTimestepLength() > couplingScheme->getTimestepLength()){
      setTimestepLength(couplingScheme->getTimestepLength());
    }
    if (getComputedTimestepPart() > couplingScheme->getComputedTimestepPart()){
      setComputedTimestepPart(couplingScheme->getComputedTimestepPart());
    }
  }
}

void CompositionalCouplingScheme:: finalize()
{
  preciceTrace("finalize()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->finalize();
  }
}

std::vector<std::string> CompositionalCouplingScheme:: getCouplingPartners() const
{
  std::vector<std::string> partners;
  std::vector<std::string> subpartners;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    subpartners = couplingScheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

void CompositionalCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver )
{
  preciceTrace("sendState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->sendState(communication, rankReceiver);
  }
}

void CompositionalCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender )
{
  preciceTrace("receiveState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->receiveState(communication, rankSender);
  }
}

std::string CompositionalCouplingScheme:: printCouplingState() const
{
  std::string state;
  std::vector<std::string> partners;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (not state.empty()){
      state += "\n";
    }
    partners = couplingScheme->getCouplingPartners();
    assertion1(partners.size() == 1, partners.size());
    state += "Coupling to ";
    state += partners[0];
    state += ":\n";
    state += couplingScheme->printCouplingState();
  }
}

void CompositionalCouplingScheme:: exportState
(
  const std::string& filenamePrefix ) const
{
  preciceTrace("exportState()");
  int enumerator = 0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    std::ostringstream stream;
    stream << filenamePrefix << "_" << enumerator;
    couplingScheme->exportState(stream.str());
    enumerator++;
  }
}

void CompositionalCouplingScheme:: importState
(
  const std::string& filenamePrefix )
{
  preciceTrace("importState()");
  int enumerator = 0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    std::ostringstream stream;
    stream << filenamePrefix << "_" << enumerator;
    couplingScheme->importState(stream.str());
    enumerator++;
  }
}

}} // namespace precice, cplscheme

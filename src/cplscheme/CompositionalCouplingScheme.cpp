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
}

void CompositionalCouplingScheme:: initialize (
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
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->addComputedTime(timeToAdd);
  }
  setTime(_couplingSchemes[0]->getTime()); // Assume all have equal time
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

std::vector<std::string> CompositionalCouplingScheme:: getCouplingPartners
(
  const std::string& accessorName ) const
{
  std::vector<std::string> partners;
  std::vector<std::string> subpartners;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    subpartners = couplingScheme->getCouplingPartners(accessorName);
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
  // TODO How to unite states?
}

void CompositionalCouplingScheme:: exportState(io::TXTWriter& writer) const
{
  preciceTrace("exportState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    // TODO Hand down filename, not writer.
    //couplingScheme->exportState(writer);
  }
}

void CompositionalCouplingScheme:: importState(io::TXTReader& reader)
{
  preciceTrace("importState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    // TODO Hand down filename, not reader.
    //couplingScheme->importState(writer);
  }
}

}} // namespace precice, cplscheme

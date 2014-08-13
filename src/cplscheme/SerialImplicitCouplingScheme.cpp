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
  SerialCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
		       secondParticipant,localParticipant,communication, dtMethod, Implicit, maxIterations)
{}

}} // namespace precice, cplscheme

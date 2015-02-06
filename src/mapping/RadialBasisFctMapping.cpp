// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "mapping/RadialBasisFctMapping.hpp"

namespace precice {
namespace mapping {

tarch::logging::Log InverseMultiquadrics:: _log("precice::mapping::InverseMultiquadrics");
tarch::logging::Log Gaussian:: _log("precice::mapping::Gaussian");
tarch::logging::Log CompactThinPlateSplinesC2:: _log("precice::mapping::CompactThinPlateSplinesC2");
tarch::logging::Log CompactPolynomialC0:: _log("precice::mapping::CompactPolynomialC0");
tarch::logging::Log CompactPolynomialC6:: _log("precice::mapping::CompactPolynomialC6");

}} // namespace precice, mapping

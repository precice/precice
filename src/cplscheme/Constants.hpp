// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_CONSTANTS
#define PRECICE_CPLSCHEME_CONSTANTS

#include <string>

namespace precice {
namespace cplscheme {
namespace constants {

const std::string& actionWriteIterationCheckpoint();

const std::string& actionReadIterationCheckpoint();

const std::string& actionWriteInitialData();

enum TimesteppingMethod
{
  FIXED_DT,
  FIRST_PARTICIPANT_SETS_DT
};


}}} // namespace precice, cplscheme, constants

#endif // PRECICE_CPLSCHEME_CONSTANTS
